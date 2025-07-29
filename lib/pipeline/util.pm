package pipeline::util;

use 5.020002;
use strict;
use warnings;
use Vcf;
use Data::Dumper;
use Carp qw(verbose confess carp croak);
require Exporter;
use List::Util qw/sum/;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use pipeline-util ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
#our %EXPORT_TAGS = ( 'all' => [ qw(
#	'TargetVcfReAnnotator'
#) ] );

#our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	LinePosLastHandleLogger
	TargetVcfReAnnotator
	WalkToTarget
	NewWalk
	WalkToNext
	AnnotateTargetRecords
	ADFilterTargetRecords
	FormatWalkTargetLineAsVcfLine
	FormatWalkTargetLineAsVcfHeader
	_formatwalkasvcflineswithfile
);

our $VERSION = "0.8.21";

# Preloaded methods go here.
# Below is stub documentation for your module. You'd better edit it!

sub TargetVcfReAnnotator{
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['targetvcf', 'vcflist']);
	my $targetvcf = Vcf->new(file=> $self -> {'targetvcf'});
	$targetvcf->parse_header();
	print $targetvcf->format_header();
	
	
	my $vcfs;
	for my $vcf (@{$self -> {'vcflist'}}){
		push(@{$vcfs},{'handle'=>Vcf->new(file=> $vcf)});
	        $vcfs -> [-1] -> {handle}->parse_header();
		$vcfs -> [-1] -> {buffer} = {'current' => [],'next' => []};
	}
	
	#iterate target file and collect data for one position across all files
	my $targetposbuffer = {'current' => [],'next' => []};
	
	my $vcfbuffers;
	#warn Dumper($targetvcf->next_data_hash()). ' ';
	while ($targetposbuffer=GetBufferedPosNext('vcf' =>\$targetvcf,'buffer' => $targetposbuffer)){
		#my $loc = $targetposbuffer -> {'current'} -> [0];
		
		for my $vcf (@{$vcfs}){
			$vcf -> {buffer} = GetBufferedPosByLoc('vcf' => \$vcf -> {handle},'buffer' => $vcf -> {buffer}, 'loc'=>  $targetposbuffer -> {'current'} -> [0]);
			#do some annotations
			#warn 'target='.
			#	$targetposbuffer -> {'current'} -> [0] -> {'CHROM'}.':'.
			#	$targetposbuffer -> {'current'} -> [0] -> {'POS'}.
			#	" vcf=".$vcf -> {'buffer'} -> {'current'} -> [0] -> {'CHROM'}.':'.
			#	$vcf -> {'buffer'} -> {'current'} -> [0] -> {'POS'}
		}
		 
	}
	#dump rest of buffer
	#insert code

}
sub NewWalk{
	#run like 	$walkdata = NewWalk('targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}}, copy_annot =>0);
	# this opens the relavant files and stores the data 

	my $self ;$self -> {'copy_annot'}=0; %{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['targetvcf', 'vcflist']);
	my $walk;
        $walk -> {'targetvcf' } -> {'file' } = $self -> {'targetvcf'};
	$walk -> {'targetvcf' } -> {'handle' } = Vcf->new(file=> $self -> {'targetvcf'});
	$walk -> {'targetvcf' } -> {'handle' } -> parse_header();
	$walk -> {'targetvcf' } -> {'buffer' } = {'current' => [],'next' => [ $walk -> {'targetvcf' } -> {handle} -> next_data_hash() ] };#GetBufferedPosNext('vcf' =>\$targetvcf,'buffer' => $targetposbuffer) or die "Empty vcf!! Cannot parse";

	#$walk -> {'vcfs' };
	for my $vcf (sort(@{$self -> {'vcflist'}})){
		push(@{$walk -> {'vcfs' }},{'file' => $vcf,'handle'=>Vcf->new(file=> $vcf)});
	        $walk -> {'vcfs' } -> [-1] -> {handle}->parse_header();
			#the folllowing block copies annotations to header if needed. 
			if($self -> {'copy_annot'}){
				for my $field (("INFO", "FORMAT")){
					for my $header (sort(keys(%{$walk -> {'vcfs' } -> [-1] -> {'handle'}->{'header'}->{$field}}))){
						if(not(defined($walk -> {'targetvcf' } -> {'handle' } -> {$field} -> {$header}))){
							#add header line here
							$walk -> {'targetvcf' } -> {'handle' } -> add_header_line(
							$walk -> {'vcfs' } -> [-1] -> {'handle'} -> {'header'} -> {$field}-> {$header});
							#append self to blabla
							#push(@{$walk -> {'targetvcf' } -> {'handle' } -> {'header_lines'}}, ${$walk -> {'targetvcf' } -> {'handle' } -> {'header'} -> {$field} -> {$header}});
						}
					}
				} 
			}
			#die Dumper($walk -> {'targetvcf' } -> {'handle' });
		$walk -> {'vcfs' } -> [-1] -> {buffer} = {'current' => [],'next' => [$walk -> {'vcfs' } -> [-1] -> {handle} -> next_data_hash()]};
	}
	
	return $walk;
}
sub WalkToNext {
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['walk']);
	my $walk = $self -> {'walk'};
	#get new position from vcf
	if(not(scalar(@{$walk -> {'targetvcf' } -> {'buffer' } -> {'next'}}) == 0)){
		$walk -> {'targetvcf' } -> {'buffer' } = GetBufferedPosNext('vcf' =>$walk -> {'targetvcf' } -> {'handle' },'buffer' => $walk -> {'targetvcf' } -> {'buffer' }) or return undef;
	}else {
		return undef;
	}
	#sync the other vcfs
	for my $vcf (@{$walk -> {'vcfs'}}){
		SyncVcfToTarget('vcf' => $vcf,'pos' => $walk -> {'targetvcf' } -> {'buffer'} -> {'current'} -> [0]);
	}
	_formatwalkasvcflineswithfile('walk' => $walk);
	if(not scalar(@{$walk -> {'targetvcf' } -> {'buffer' } -> {'next'}})){
		return undef;
	}
	return 1;
}

sub AnnotateTargetRecords {
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['walk']);
	my $walk = $self -> {'walk'};
	my $annotated;
	#die Dumper($self -> {'record'}).  "  ";
	for my $target ($walk -> {'targetvcf' }){
		for my $targetrecord (@{$target -> {'buffer' } -> {'current'}}){
			for my $vcf (@{$walk -> {'vcfs' }}){
				for my $record (@{$vcf -> {'buffer' } -> {'current'}}){
					#die Dumper($self -> {'record'}).  "  ";
					next if(not(VariantIsEq('loc1' => $targetrecord, 'loc2' => $record)));
					
					#annotate
					#$walk -> {'targetvcf' }  -> {'buffer' } -> {'current'} = 
					$targetrecord = AnnotateVariant('targetvcfhandle' => $walk -> {'targetvcf' } -> {'handle'},'targetrecord' => $targetrecord, 'record' => $record);
					$annotated -> { $vcf -> {'file'} }++; 
				}
			}
			#Dumper($record);
			#next if not(scalar(keys(%{$record})));
			#$out .= $vcf -> {'file' } . ":\tcurr:\t" .
			#	$vcf -> {'handle'} -> format_line($record);
		}
	}
	return $walk;
}

sub ADFilterTargetRecords {
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['walk','deltafrequency','deltacount']);
	my $walk = $self -> {'walk'};
	my $annotated;
	#die Dumper($self -> {'record'}).  "  ";
	for my $target ($walk -> {'targetvcf' }){
		for my $targetrecord (@{$target -> {'buffer' } -> {'current'}}){
			for my $vcf (@{$walk -> {'vcfs' }}){
				for my $record (@{$vcf -> {'buffer' } -> {'current'}}){
					#die Dumper($self -> {'record'}).  "  ";
					#warn Dumper({'loc1' => $targetrecord, 'loc2' => $record,"current" => @{$vcf -> {'buffer' } -> {'current'}},"next" => @{$vcf -> {'buffer' } -> {'current'}}}) if(not(VariantIsEq('loc1' => $targetrecord, 'loc2' => $record)));
					#next if(not(VariantIsEq('loc1' => $targetrecord, 'loc2' => $record)));
					
					#annotate
					#$walk -> {'targetvcf' }  -> {'buffer' } -> {'current'} = 
					$targetrecord = ADFilterVariant('targetvcfhandle' => $walk -> {'targetvcf' } -> {'handle'},'targetrecord' => $targetrecord, 'record' => $record,'deltafrequency'=> $self -> {'deltafrequency'},'deltacount'=> $self -> {'deltacount'});
					$annotated -> { $vcf -> {'file'} }++; 
				}
			}
			#Dumper($record);
			#next if not(scalar(keys(%{$record})));
			#$out .= $vcf -> {'file' } . ":\tcurr:\t" .
			#	$vcf -> {'handle'} -> format_line($record);
		}
	}
	return $walk;
}

sub AnnotateVariant {
	my $self ;%{$self}= @_;
	SelfRequire(%{$self},'req'=> ['targetrecord','targetvcfhandle','record']);
	#die Dumper($self -> {'record'}).  "  ";
	for my $sample (keys(%{$self -> {'targetrecord'} -> {'gtypes'}})){
		
		if(defined($self -> {'record'} -> {'gtypes'} -> {$sample})){
			#die Dumper($self -> {'record'})."object dump" ;
			for my $infofield (sort(keys(%{$self -> {'record'} ->  {'INFO'}}))){
				#second part needs to be either . or .(,.)* so to make it fast substr($text,0,1) eq '.'
				if(not(defined $self -> {'targetrecord'} -> {'INFO'} -> {$infofield}) || substr($self -> {'targetrecord'} -> {'INFO'} -> {$infofield},0,1) eq '.'){
					$self -> {'targetrecord'} -> {'INFO'} -> {$infofield} = $self -> {'record'} -> {'INFO'} -> {$infofield};
				}
				
			}
			#warn Dumper(sort(keys(%{$self -> {'record'} ->  {'gtypes'} -> {$sample} })))." ################################################################################################# ";
			for my $formatfield (sort(keys(%{$self -> {'record'} ->  {'gtypes'} -> {$sample} }))){

				#This might contain skips if for example the genotype is already present
				if(defined($formatfield) && $formatfield eq 'GT' && 
					($self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} eq '.' ||
					 $self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} eq './.') ){
					$self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} = $self -> {'record'} ->  {'gtypes'} -> {$sample} -> {$formatfield};
				}
				#second part needs to be either . or .(,.)* so to make it fast substr($text,0,1) eq '.'
				if( (not( defined($self -> {'targetrecord'} ->  {'gtypes'} -> {$sample} -> {$formatfield})) || 
					 substr($self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield},0,1) eq '.') && 
					 defined($self -> {'record'} -> {'gtypes'} -> {$sample} -> {$formatfield}) &&
					 substr($self -> {'record'} -> {'gtypes'} -> {$sample} -> {$formatfield},0,1) ne '.'){
					
					$self -> {'targetvcfhandle'} -> add_format_field($self -> {'targetrecord'},$formatfield);
					$self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} = $self -> {'record'} -> {'gtypes'} -> {$sample} -> {$formatfield};
					#$vcf->add_format_field($x,'FOO');
				}
			}
		}
	}
	return $self -> {'targetrecord'};
}

sub ADFilterVariant {
	#filter by ad count and frequency
	# code is messy but should run efficient

	my $self;%{$self}= @_;
	SelfRequire(%{$self},'req'=> ['targetrecord','targetvcfhandle','record']);
	#die Dumper($self -> {'record'}).  "  ";
	my $failsNoAdFilter = 0;
	my $failsCountFilter = 0;
	my $failsFreqFilter = 0;	
	#warn Dumper($self -> {'targetrecord'},$self -> {'record'})." " ;
	#process AD by target samples in vcf
	for my $targetsample (keys(%{$self -> {'targetrecord'} -> {'gtypes'}})){
		#warn "The ad is seen as:"._isValidAD($self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'});
		if(not(defined($self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}))){
					#warn "No AD target record";
					$failsNoAdFilter++;
		}elsif(not(_isValidAD($self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}))){
					#warn "No AD target record";
					$failsNoAdFilter++;
		}else{
			#warn "Chromposiseq".ChromPosIsEq('loc1' => $self -> {'targetrecord'}, 'loc2' => $self -> {'record'});
			if(ChromPosIsEq('loc1' => $self -> {'targetrecord'}, 'loc2' => $self -> {'record'}) &&
				$self -> {'targetrecord'} -> {'REF'} eq $self -> {'record'} -> {'REF'}){
				
				#warn " ## ### Comparison ".Dumper($self -> {'targetrecord'}  -> {'CHROM'},$self -> {'targetrecord'}  -> {'POS'})."vs". Dumper($self -> {'record'} -> {'CHROM'}, $self -> {'record'} -> {'POS'});
				
				#filtersample: all the samples to filter against AD count/frequency + observed frequency seen in samples
				for my $filtersample (keys(%{$self -> {'record'} -> {'gtypes'}})){
					if($targetsample eq $filtersample){
						# case to protect filtering against self
						#next
					}else{
					
						#vcf oddities on empty targetsample fields
						if(not(defined($self -> {'record'} -> {'gtypes'} -> {$filtersample} -> {'AD'}))){
							if(FailsSimpleCountFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'mincount' => $self -> {'deltacount'})){
								$failsCountFilter++;
								#warn " ## ### Comparison ";
							};
							if(FailsSimpleFreqFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'minfreq' => $self -> {'deltafrequency'})){
								$failsFreqFilter++;
								#warn " ## ### Comparison ";
							};
						}elsif(not(_isValidAD($self -> {'record'} -> {'gtypes'} -> {$filtersample} -> {'AD'}))){
						#warn " ## ### Comparison ";
							if(FailsSimpleCountFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'mincount' => $self -> {'deltacount'})){
								$failsCountFilter++;
								#warn " ## ### Comparison ";
							};
							if(FailsSimpleFreqFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'minfreq' => $self -> {'deltafrequency'})){
								$failsFreqFilter++;
								#warn " ## ### Comparison ";
							};
						}else{
							#warn " ## ### Comparison ";
							#both should be present
							#warn "both should be present".ChromPosIsEq('loc1' => $self -> {'targetrecord'}, 'loc2' => $self -> {'record'});
							#there 3 if statments are to speed things up
							if(FailsSimpleCountFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'mincount' => $self -> {'deltacount'})){
								$failsCountFilter++;
								#warn " ## ### Comparison ";
							};
							if(FailsSimpleFreqFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'minfreq' => $self -> {'deltafrequency'})){
								$failsFreqFilter++;
								#warn " ## ### Comparison ";
							};
							if($failsCountFilter == 0 || $failsFreqFilter == 0){
								#warn " ## ### Comparison ";
								#warn "## DELTAFILTER".ChromPosIsEq('loc1' => $self -> {'targetrecord'}, 'loc2' => $self -> {'record'});
								#$passNoAdFilter++;#should be filtered???
								#defined($self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}) or die "Error here".Dumper($self -> {'targetrecord'}). " ";
								#defined($self -> {'record'} -> {'gtypes'} -> {$filtersample} -> {'AD'}) or die "Error here".Dumper($self -> {'record'},$self -> {'targetrecord'}). " ";
								my @targetAD=split(',',$self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}) or die "Error here".Dumper($self -> {'targetrecord'}). " ";
								my @targetALT= @{$self -> {'targetrecord'} -> {'ALT'}};
								my @filterAD=split(',',$self -> {'record'} -> {'gtypes'} -> {$filtersample} -> {'AD'}) or die "Error here".Dumper($self -> {'record'}). " ";
								my @filterALT= @{$self -> {'record'} -> {'ALT'}};

								#warn "No AD and filter AD record";

								#next if($self -> {'record'} -> {'gtypes'} -> {$filtersample} -> {'AD'} eq '.' || 
								#	$self -> {'record'} -> {'gtypes'} -> {$filtersample} -> {'AD'} eq ".".",." x scalar(@{$self -> {'record'} -> {'ALT'}}));

								my $targetindex = 0;
								while ($targetindex < scalar(@targetALT)){
									#warn "loop1";
									
									my $filterindex = 0;
									#should be cleaned with my $filterindex = GetFilterIndex('alt'=>$targetAlt[$targetindex],'filteralt'=> \@filterALT); 
									while ($filterindex < scalar(@filterALT)){
										#warn "loop2";
										
										#warn "#!!Comparison ".$targetALT[$targetindex ]."vs". $filterALT[$filterindex]. 
										#	" count ".$targetAD[$targetindex + 1 ].">=".($self -> {'deltacount'} + $filterAD[$filterindex + 1 ]);
										#warn "#!!Comparison ".$targetALT[$targetindex ]."vs". $filterALT[$filterindex]. 
										#	" freq ".$targetAD[$targetindex + 1 ]/sum(@targetAD).">=". ($self -> {'deltafrequency'} + $filterAD[$filterindex + 1]/sum(@filterAD));
										#warn "#!!Comparison ".$filterAD[$filterindex + 1].sum(@filterAD).($self -> {'deltafrequency'} + $filterAD[$filterindex + 1]/sum(@filterAD));
										#filter for count
										if($targetALT[$targetindex] eq $filterALT[$filterindex] && $targetAD[$targetindex  + 1] < ($self -> {'deltacount'} + $filterAD[$filterindex + 1]) ){
											$failsCountFilter++;
										}
										#filter for freq
										if($targetALT[$targetindex] eq $filterALT[$filterindex] && 
											($targetAD[$targetindex + 1]/sum(@targetAD)) < ($self -> {'deltafrequency'} + $filterAD[$filterindex + 1]/sum(@filterAD)) ){						
											$failsFreqFilter++;
										}
										$filterindex++;
									}
									$targetindex++;
								}
								
							}
							
						}
					}
				}
			}else{
				#warn " Comparison chromposref not matched ".Dumper($self -> {'targetrecord'}  -> {'CHROM'},$self -> {'targetrecord'}  -> {'POS'})."vs". Dumper($self -> {'record'} -> {'CHROM'}, $self -> {'record'} -> {'POS'});
				if(FailsSimpleCountFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'mincount' => $self -> {'deltacount'})){
					$failsCountFilter++;
				};
				if(FailsSimpleFreqFilter('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'minfreq' => $self -> {'deltafrequency'})){
					$failsFreqFilter++;
				};
			}
		}
	}
	if($self -> {'targetrecord'} -> {'FILTER'} -> [0] eq '.'){
		$self -> {'targetrecord'} -> {'FILTER'} -> [0] = 'PASS';
	}
	#if($passNoAdFilter == 0 && index(join(";",@{$self -> {'targetrecord'} -> {'FILTER'}}),"NoAD") == -1){
#
#		if($self -> {'targetrecord'} -> {'FILTER'} -> [0] eq 'PASS' || $self -> {'targetrecord'} -> {'FILTER'} -> [0] eq '.'){
#			@{$self -> {'targetrecord'} -> {'FILTER'}} = ("NoAD");	
#		}else{
#			push(@{$self -> {'targetrecord'} -> {'FILTER'}},"NoAD");
#		}
#	}
	if($failsCountFilter > 0 && index(join(";",@{$self -> {'targetrecord'} -> {'FILTER'}}),'ADDeltaCountlt'.$self -> {'deltacount'}) == -1){
		if($self -> {'targetrecord'} -> {'FILTER'} -> [0] eq 'PASS' || $self -> {'targetrecord'} -> {'FILTER'} -> [0] eq '.'){							
			@{$self -> {'targetrecord'} -> {'FILTER'}} = ('ADDeltaCountlt'.$self -> {'deltacount'});	
		}else{
			push(@{$self -> {'targetrecord'} -> {'FILTER'}},'ADDeltaCountlt'.$self -> {'deltacount'});
		}
	}
	if($failsFreqFilter > 0 && index(join(";",@{$self -> {'targetrecord'} -> {'FILTER'}}),"ADDeltaFrequencylt".$self -> {'deltafrequency'}) == -1){
		if($self -> {'targetrecord'} -> {'FILTER'} -> [0] eq 'PASS' || $self -> {'targetrecord'} -> {'FILTER'} -> [0] eq '.'){
			@{$self -> {'targetrecord'} -> {'FILTER'}} = ("ADDeltaFrequencylt".$self -> {'deltafrequency'});	
		}else{
			push(@{$self -> {'targetrecord'} -> {'FILTER'}},"ADDeltaFrequencylt".$self -> {'deltafrequency'});
		}
	}
	return $self -> {'targetrecord'};
}
#			#die Dumper($self -> {'record'})."object dump" ;
#			for my $infofield (sort(keys(%{$self -> {'record'} ->  {'INFO'}}))){
#				#second part needs to be either . or .(,.)* so to make it fast substr($text,0,1) eq '.'
#				if(not(defined $self -> {'targetrecord'} -> {'INFO'} -> {$infofield}) || substr($self -> {'targetrecord'} -> {'INFO'} -> {$infofield},0,1) eq '.'){
#					$self -> {'targetrecord'} -> {'INFO'} -> {$infofield} = $self -> {'record'} -> {'INFO'} -> {$infofield};
#				}
#				
#			}
#			#warn Dumper(sort(keys(%{$self -> {'record'} ->  {'gtypes'} -> {$sample} })))." ################################################################################################# ";
#			for my $formatfield (sort(keys(%{$self -> {'record'} ->  {'gtypes'} -> {$sample} }))){
#
#				#This might contain skips if for example the genotype is already present
#				if(defined($formatfield) && $formatfield eq 'GT' && 
#					($self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} eq '.' ||
#					 $self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} eq './.') ){
#					$self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} = $self -> {'record'} ->  {'gtypes'} -> {$sample} -> {$formatfield};
#				}
#				#second part needs to be either . or .(,.)* so to make it fast substr($text,0,1) eq '.'
#				if( (not( defined($self -> {'targetrecord'} ->  {'gtypes'} -> {$sample} -> {$formatfield})) || 
#					 substr($self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield},0,1) eq '.') && 
#					 defined($self -> {'record'} -> {'gtypes'} -> {$sample} -> {$formatfield}) &&
#					 substr($self -> {'record'} -> {'gtypes'} -> {$sample} -> {$formatfield},0,1) ne '.'){
#					
#					$self -> {'targetvcfhandle'} -> add_format_field($self -> {'targetrecord'},$formatfield);
#					$self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} = $self -> {'record'} -> {'gtypes'} -> {$sample} -> {$formatfield};
#					#$vcf->add_format_field($x,'FOO');
#				}
#			}

sub _isValidAD {
	my $AD;
	$AD = $_[0] or confess "Needs string AD field input.";
	my @DepthByAlle= split(',',$AD);
	#warn Dumper(sum( map { if($_ eq '.'){$_ = 0}   } @DepthByAlle))." ";
	if(not(defined($AD))){
		return 0;
	}elsif($AD eq '.'|| $AD eq ".".",." x scalar(@DepthByAlle)){
		return 0;
	}else{
		map { if($_ eq "."){$_ = 0} }(@DepthByAlle);
		if(sum(@DepthByAlle) <= 0){
			return 0;
		}
	}
	return 1;
}
sub FailsSimpleCountFilter {
	#('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'mincount' => ,'minfreq' => $self -> {'deltacount'}, $self -> {'deltafrequency'})
	my $self;%{$self}= @_;
	SelfRequire(%{$self},'req'=> ['ad','mincount']);
	my @targetAD=split(',',$self -> {'ad'});
	#my @targetALT= @{$self -> {'targetrecord'} -> {'ALT'}};
	
	#warn "No AD Filter record";
	#should skip the ref count and AD is assumed to be per possible allele aka REF + ALT fields and not per Gentype...
	my $targetindex = 1;#skips the ref field
	while ($targetindex < (scalar(@targetAD))){
		#warn "Comparison ". 
		#	" count ".$targetAD[$targetindex ].">=".($self -> {'mincount'});
		if($targetAD[$targetindex ] < ($self -> {'mincount'})){
			return 1;
		};
		#if($targetAD[$targetindex + 1 ]/sum(@targetAD) >= ($self -> {'deltafrequency'})){
		#	$passFreqFilter++;
		#};
		$targetindex++;
	}
	return 0;
	
}
sub FailsSimpleFreqFilter {
	#('ad'=> $self -> {'targetrecord'} -> {'gtypes'} -> {$targetsample} -> {'AD'}, 'mincount' => ,'minfreq' => $self -> {'deltacount'}, $self -> {'deltafrequency'})
	my $self;%{$self}= @_;
	SelfRequire(%{$self},'req'=> ['ad','minfreq']);
	my @targetAD=split(',',$self -> {'ad'});
	#my @targetALT= @{$self -> {'targetrecord'} -> {'ALT'}};
	
	#warn "No AD Filter record";
	#should skip the ref count and AD is assumed to be per possible allele aka REF + ALT fields and not per Gentype...
	my $targetindex = 1;#skips the ref field
	while ($targetindex < (scalar(@targetAD))){
		#if($targetAD[$targetindex ] <= ($self -> {'deltacount'})){
		#	return 0;
		#};
		#warn "Comparison ".
		#	" freq ".$targetAD[$targetindex ]/sum(@targetAD).">=".$self -> {'minfreq'};
		if($targetAD[$targetindex ]/sum(@targetAD) < ($self -> {'minfreq'})){
			return 1;
		};
		$targetindex++;
	}
	return 0;
}
sub _formatwalkasvcflineswithfile {
	#Quick dump of walk data for debugging;
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['walk']);
	my $walk = $self -> {'walk'};
	#$vcf->format_line($x);
	my $out;	
	$out .= '<--------'."\n"; 
	for my $vcf ($walk -> {'targetvcf' },@{$walk -> {'vcfs' }}){
		for my $record (@{$vcf -> {'buffer' } -> {'current'}}){
			#die Dumper($record);
			next if not(scalar(keys(%{$record})));
			$out .= $vcf -> {'file' } . ":\tcurr:\t" .
				$vcf -> {'handle' } -> format_line($record);
		}
		for my $record (@{$vcf -> {'buffer' } -> {'next'}}){
			next if not(scalar(keys(%{$record})));
			$out .= $vcf -> {'file' } . ":\tnext:\t" .
				$vcf -> {'handle' } -> format_line($record);
		}
	}
	$out .= '-------->'."\n ";
	
	return $out;
	
}

sub FormatWalkTargetLineAsVcfLine {
	#Output Target vcf data;
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['walk']);
	my $walk = $self -> {'walk'};
	#$vcf->format_line($x);
	my $out;	
	#$out .= '<-----'."\n"; 
	for my $vcf ($walk -> {'targetvcf' }){
		for my $record (@{$vcf -> {'buffer' } -> {'current'}}){
			#Dumper($record);
			next if not(scalar(keys(%{$record})));
			$out .= $vcf -> {'handle' } -> format_line($record);
		}
	}
	#$out .= '----->'."\n";
	
	return $out;
}
sub FormatWalkTargetLineAsVcfHeader {
	#Output Target vcf header;
	
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['walk']);
	my $walk = $self -> {'walk'};
		my $out;
	
	for my $header (@{$self -> {'headerlines'} }){
		#warn Dumper($header)." ";
		$walk -> {'targetvcf' } -> {'handle' } -> add_header_line(\%{$header});
	}
	$out .= $walk -> {'targetvcf' } -> {'handle' } -> format_header();
	#$out .= 
	
	return $out;
	
}

sub SyncVcfToTarget{
	my $self ;%{$self}= @_ ;
	
	#@{$self -> {'vcflist'}}=($self -> {'targetvcf'});
#	SelfRequire(%{$self},'req'=> ['walk']);
	SelfRequire(%{$self},'req'=> ['vcf','pos']);	
	#my $walk = $self -> {'walk'};
	my $continue = 1;
	return if(not(defined($self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0])));
	if(not(ChromPosIsNext('curr' => GetLoc($self -> {'pos'}), 'next' => GetLoc($self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0]), 'vcf' => $self -> {'vcf' } -> {'handle'}))){
		while (($continue == 1) && ($self -> {'vcf' } -> {'buffer'} = GetBufferedPosNext('vcf' =>$self -> {'vcf' } -> {'handle'},'buffer' => $self -> {'vcf' } -> {'buffer'}))){
			
			#Goes until equal position or greater then position; #not(ChromPosIsNext('loc1' => {'CHROM'=>4,'POS'=>55593536}, 'loc2' => $targetposbuffer -> {'current'} -> [0], 'vcf' => $targetvcf))
			
			#Goes till
			if(not(defined($self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0])) || ChromPosIsNext('curr' => GetLoc($self -> {'pos'}), 'next' => GetLoc($self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0]), 'vcf' => $self -> {'vcf'} -> {'handle'})){
				
				#warn LocGetChromPosAsString(loc => GetLoc($self -> {'vcf'} -> {'buffer'} -> {'current'} -> [0]))." ";
				
				$continue = 0;
				
			}
			
		}
	}
	#if( scalar(keys %{$self -> {'vcf'} -> {'buffer'} -> {'current'} -> [0]})){
	#	warn "Current ".Dumper(GetLoc(%{$self -> {'vcf'} -> {'buffer'} -> {'current'} -> [0]})). " ";
	#}
	#if( scalar(keys %{$self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0]}) ){
	#	warn "Next ".Dumper(GetLoc(%{$self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0]})). " ";
	#}
	#warn "Sync done ";
}

sub WalkToTarget{
	my $self ;%{$self}= @_ ;
	@{$self -> {'vcflist'}}=($self -> {'targetvcf'});
	SelfRequire(%{$self},'req'=> ['targetvcf', 'vcflist','pos']);
	my $targetvcf = Vcf->new(file=> $self -> {'targetvcf'});
	$targetvcf->parse_header();
	print $targetvcf->format_header();

	
	my $vcfs;
	for my $vcf (@{$self -> {'vcflist'}}){
		push(@{$vcfs},{'handle'=>Vcf->new(file=> $vcf)});
	        $vcfs -> [-1] -> {handle}->parse_header();
		$vcfs -> [-1] -> {buffer} = {'current' => [],'next' => []};
	}
	
	#iterate target file and collect data for one position across all files
	my $targetposbuffer = {'current' => [],'next' => []};
	
	my $vcfbuffers;
	#warn Dumper($targetvcf->next_data_hash()). ' ';
	my $continue = 1;
	while ($continue == 1 && ($targetposbuffer=GetBufferedPosNext('vcf' =>\$targetvcf,'buffer' => $targetposbuffer))){
		
		#Goes until equal position or greater then position; #not(ChromPosIsNext('loc1' => {'CHROM'=>4,'POS'=>55593536}, 'loc2' => $targetposbuffer -> {'current'} -> [0], 'vcf' => $targetvcf))
		
		#Goes till
		if(ChromPosIsNextOrEqual('curr' => $self -> {'pos'}, 'next' => $targetposbuffer -> {'next'} -> [0], 'vcf' => $targetvcf)){
			
			#warn Dumper({'CHROM' => $targetposbuffer -> {'current'} -> [0] -> {'CHROM'},'POS' => $targetposbuffer -> {'current'} -> [0] -> {'POS'}}) . "  ";
			
			$continue = 0;
			
		}
		
	}
	#warn Dumper(GetLoc(%{$targetposbuffer -> {'current'} -> [0]}));
	#warn "done ";
}
sub SelfRequire {
        #calls like SelfRequire(%{$self}, 'req'=> ['vcf', 'fields']);
        my $self;#carp "## ";
        %{$self}=@_ or  confess "input not parseable as hash:".Dumper(@_). ' ';
        #die Dumper($self)." ";
        my $reqs = $self -> {'req'} or die "no requirement array as input";
        for my $tag (@{$reqs}){
                defined($self -> {$tag}) or confess "no requirement '$tag' as input";
        }
}
sub GetBufferedPosNext {
	my $self; %{$self}=@_;
	SelfRequire(%{$self}, 'req'=> ['vcf', 'buffer']);
	#warn Dumper($self -> {'vcf'}).' ';
	my $targetvcf = $self -> {'vcf'} ;#or confess("Error.");
	         #warn Dumper($targetvcf->next_data_hash()). ' ';
	
	$self -> {'buffer'} = PosBufferBufferCleanUp('buffer' => $self -> {'buffer'});
	if(scalar(@{$self -> {'buffer'} -> {'next'}})){
		$self -> {'buffer'} -> {'current'} = $self -> {'buffer'} -> {'next'};
		$self -> {'buffer'} -> {'next'} = [];
	}
	#what works : defined($self -> {'buffer'} -> {'current'}) or defined($self -> {'buffer'} -> {'current'} -> [0] ? )
	#warn '## Buffer ##'. Dumper($self -> {'buffer'})." ";
	my $x;
	#warn "next amount ".scalar(@{$self -> {'buffer'} -> {'next'}})." ";
	while( scalar(@{$self -> {'buffer'} -> {'next'}}) == 0 and my $x=$targetvcf -> next_data_hash()){
	#while( scalar(@{$self -> {'buffer'} -> {'next'}}) == 0 and defined($x)){
	#	$x=$targetvcf -> next_data_hash();
		#warn '## VC ## '. Dumper($x). ' ';
		push(@{$self -> {'buffer'} -> {'next'}},$x);
		$self -> {'buffer'} = PosBufferBufferCleanUp('buffer' => $self -> {'buffer'});

	}


		#croak(Dumper($self -> {'buffer'})) if (@{$self -> {'buffer'} -> {'current'}} >1);
		#carp "x=".Dumper($x)."vcf=".Dumper($targetvcf).".result## ".scalar(@{$self -> {'buffer'} -> {'current'}})." ";
	if(scalar(@{$self -> {'buffer'} -> {'current'}}) > 0 ){
		return $self -> {'buffer'} ;
	}else{
		return undef;
	}
}
sub PosBufferBufferCleanUp {
	my $self; %{$self}=@_;
	#checks if the value should be added to current next or return undef??
	SelfRequire(%{$self}, 'req'=> [ 'buffer']);
	if(scalar(@{$self -> {'buffer'} -> {'current'}}) == 0){
		#init and add first val of array
		$self -> {'buffer'} -> {'current'} = $self -> {'buffer'} -> {'next'} or confess Dumper($self -> {'buffer'} -> {'next'});
		$self -> {'buffer'} -> {'next'} = [];
		#push(@{$self -> {'buffer'} -> {'current'}},$x);
		#warn 'Init';
	}elsif(scalar(@{$self -> {'buffer'} -> {'current'}}) && defined($self -> {'buffer'} -> {'next'} -> [0] -> {'POS'}) && 
		$self -> {'buffer'} -> {'next'} -> [0] -> {'POS'} eq $self -> {'buffer'} -> {'current'} -> [0] ->  {'POS'} && 
		$self -> {'buffer'} -> {'next'} -> [0] -> {'CHROM'} eq $self -> {'buffer'} -> {'current'} -> [0] ->  {'CHROM'}){
		#warn 'Grow';
		#grow array
		push(@{$self -> {'buffer'} -> {'current'}},@{$self -> {'buffer'} -> {'next'}});
		$self -> {'buffer'} -> {'next'} = [];
	}else{
		#chrom / pos not the same so dump as warning
		#croak(Dumper($self -> {'buffer'}));
		#warn 'Next';
		#assign to next
		#$self -> {'buffer'} -> {'next'} = [$x];
		#process next
	}
	return $self -> {'buffer'};
}
#sub GetBufferedPosByLoc {
#	my $self; %{$self}=@_;
#	SelfRequire(%{$self}, 'req'=> ['vcf', 'buffer','loc']);
#	
#	#$self->{'buffer'}=GetBufferedPosNext('vcf'=>$self->{'vcf'},'buffer'=>$self->{'buffer'});
#	warn LocGetChromPosAsString('loc'=>  $self -> {'loc'})." ";
#	#scalar(@{$self -> {'buffer'} -> {'next'}}) == 0)
#	my $count = 0;
#	while( (ChromPosIsNext('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0],'vcf'=>$self->{'vcf'}) || 
#		not(ChromPosIsEq('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0]))) && 
#		scalar($self->{'buffer'}->{'current'}) > 0){
#		$count++;
#		#croak "test";
#		$self->{'buffer'}=GetBufferedPosNext('vcf'=>$self->{'vcf'},'buffer'=>$self->{'buffer'});
#		last if(not(defined($self->{'buffer'} -> {current}) && ref($self->{'buffer'} -> {current}) eq 'ARRAY' ) );		
#		warn join (",\n",($count, LocGetChromPosAsString('loc'=>$self -> {'buffer'} -> {'current'} -> [0]) , LocGetChromPosAsString('loc'=>$self -> {'loc'})))." ";		
#		die "problems here" if($count > 10);
#	}
#	
#	if(not(defined($self->{'buffer'} -> {current}) && ref($self->{'buffer'} -> {current}) eq 'ARRAY' ) ){
#		return undef;
#	}elsif(ChromPosIsEq('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0])){
#		#return if equal
#		return $self -> {'buffer'};
#	}elsif(ChromPosIsNext('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0],'vcf'=>$self->{'vcf'})){
#		#return if higher then pos #note the same behavior as above, im not sure what to do here
#		return $self -> {'buffer'};
#	}else{
#		return undef;
#	}
#}
#sub GetBufferedPosByLoc {
#	
#}
sub LinePosLastHandleLogger {
	my $self; %{$self}=@_;
	SelfRequire(%{$self}, 'req'=> ['lineinterval', 'timeinterval']);
	if($. % $self -> {'modulo'} == 0){
		warn localtime(time())."\tINFO\t$0:\tIterating records currently at line $.\n" if(time % $self -> {timeinterval} == 0 );
	}
}

sub GetContigsAsArray {
	my $self; %{$self} = @_;
	SelfRequire(%{$self},'req'=> ['vcf']);
	my @contigs;
	#confess "ERROR Undefinded" .Dumper($self -> {'vcf'});
	my $vcf = $self -> {'vcf'} or confess "ERROR Undefinded" .Dumper($self -> {'vcf'});
	#warn Dumper($vcf -> {'header_lines'})." ";
	for my $hdata (@{$vcf -> {'header_lines'}}){
		if($hdata -> {'key'} eq 'contig'){
			push(@contigs, $hdata ->{'ID'});
		}
	}
	#die Dumper(@contigs);
	#die Dumper(${$self -> {vcf}} -> get_header_line('key'=>'contig'))." ";
	return @contigs;
}

sub ContigIsNext {
	#checks if 
	my $self ;%{$self} = @_;
	SelfRequire(%{$self},'req'=> ['curr','next','vcf']);
	my $contigcur=0;
	my $contigidxcurr=0;my $contigidxnext=0;
	#warn Dumper($self -> { 'vcf'}). ' ' if($. > 300);
	defined($self ->  { 'curr'} -> {'CHROM'}) or confess("No requirement curr CHROM as input" . ' ');
	defined($self ->  { 'next'} -> {'CHROM'}) or confess("No requirement next CHROM as input" . ' ');
	for my $contig (GetContigsAsArray('vcf' => $self -> { 'vcf'})){
		$contigidxcurr=$contigcur if(defined $self -> { 'curr'} -> {'CHROM'} && $self -> { 'curr'} -> {'CHROM'} eq $contig);
		$contigidxnext=$contigcur if(defined $self -> { 'next'} -> {'CHROM'} && $self -> { 'next'} -> {'CHROM'} eq $contig);
		$contigcur++;
	}
	#warn "############ $contigidxcurr > $contigidxnext";
	return 1 if($contigidxcurr < $contigidxnext );

	return 0;#or should i die here?
}

sub ChromPosIsEq {
	my $self;%{$self} = @_;
	#warn "loc1\t".LocGetChromPosAsString('loc' => GetLoc($self -> {'loc1'}))."\tloc2\t". LocGetChromPosAsString('loc' => GetLoc($self -> {'loc2'})). " ";
	if(defined($self -> { 'loc2'}) &&
		defined($self -> { 'loc2'} -> {'CHROM'}) &&
		$self -> { 'loc1'} -> {'CHROM'} eq $self -> { 'loc2'} -> {'CHROM'} && 
		$self -> { 'loc1'} -> {'POS'} eq $self -> { 'loc2'} -> {'POS'}){
		#warn "###\n######IsEq";
		return 1;
	}else{
		return 0;
	}
}

sub VariantIsEq {
	my $self;%{$self} = @_;
	#aka chrom pos ref alt is equal
	if(ChromPosIsEq(%{$self}) && $self -> { 'loc1'} -> {'REF'} eq $self -> { 'loc2'} -> {'REF'} && scalar@{$self -> { 'loc1'} -> {'ALT'}} == scalar @{$self -> { 'loc2'} -> {'ALT'}}){
		my $alteq = 1;
		my $altidx = 0;
		while($alteq == 1 && $altidx < scalar(@{$self -> { 'loc1'} -> {'ALT'}}) ){
			#Check if alts are equal
			$alteq = 0 if $self -> { 'loc1'} -> {'ALT'} -> [$altidx] ne $self -> { 'loc2'} -> {'ALT'} -> [$altidx];
			#if($annotated == 1){
			#	return 1;
			#}
			$altidx++;
		}
		return $alteq;
	}else {
		return 0;
	}
}

sub ChromPosIsNext  {
	#checks if loc1 > loc2
	my $self; %{$self} = @_;
	SelfRequire(%{$self},'req'=> ['curr','next','vcf']);
	#carp "loc1".Dumper($self -> { 'loc1'})."loc2".Dumper($self -> { 'loc2'});
	if((not(defined( $self -> { 'next'} -> {'CHROM'}))|| not(defined( $self -> { 'curr'} -> {'CHROM'})))){
		return 1;
	}elsif((($self -> { 'curr'} -> {'CHROM'} eq $self -> { 'next'} -> {'CHROM'} && 
			$self -> { 'curr'} -> {'POS'} < $self -> { 'next'} -> {'POS'})) || 
		(defined( $self -> { 'next'}) && defined( $self -> { 'curr'} -> {'CHROM'}) && $self -> { 'curr'} -> {'CHROM'} ne $self -> { 'next'} -> {'CHROM'} && ContigIsNext(%$self)) ){
		#warn '###########ChromPosIsNext::ret='.1;
		# next position is next to current posiotion that the files are being synced to so say ok and stop reading this file...
		return 1;
	}else{
		# read on my friend
		return 0;
	}
}

sub ChromPosIsNextOrEqual  {
	my $self; %{$self} = @_;
	SelfRequire(%{$self},'req'=> ['curr','next','vcf']);
	#carp "loc1".Dumper($self -> { 'loc1'})."loc2".Dumper($self -> { 'loc2'});
	if(ChromPosIsNext(%{$self}) || ChromPosIsEq(%{$self})){
		#warn 'ChromPosIsNext::ret='.1;
		return 1;
	}else{
		return 0;
	}
}
sub LocGetChromPosAsString{
	my $self; %{$self} = @_;
	return $self -> {'loc'} -> {'CHROM'}.":".$self -> {'loc'} -> {'POS'};
}
#sub GetLoc {
#	my $self;
#       %{$self}=@_ or  confess "input not parseable as hash:".Dumper(@_). ' ';
#
#	return {'CHROM' => $self -> {'CHROM'},'POS' => $self -> {'POS'}};
#}
sub GetLoc {
	my $self; $self = shift @_ or  confess "input not parseable as hash:".Dumper(@_). ' ';
	my $loc;%{$loc} =  ('CHROM' => $self -> {'CHROM'},
		'POS' => $self -> {'POS'});
	return $loc;
}
1;
__END__

=head1 NAME

pipeline::util - Perl cmdline tools for hts pipeline

=head1 SYNOPSIS

execute scripts from bash

=head1 DESCRIPTION

read readme on github

=head2 EXPORT

None by default.

=head1 SEE ALSO

	https://github.com/mmterpstra/pipeline-util

=head1 AUTHOR

m.m.terpstra

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2018 by m.m.terpstra

LGPL licence

=cut

=head2 TargetVcfReAnnotator
    About   : function which lifts over annotations which are lost by merges
    Usage   : TargetVcfReAnnotator('targetvcf'=>'test.vcf','vcflist'=>['a.vcf,b.vcf']);               # (from a script)
    Args    : File names of A the target vcf to annotate and B the list of vcfs to resque annotations from.
=cut


=head2 SelfRequire
    About   : minimal validation tool of subroutine input
    Usage   : SelfRequire(%{$self}, 'req'=> ['vcf', 'fields']); assumes input of self like a hash and checks it the fields are present in this hash
    Args    : 'req' = Fields to present in self hash
=cut

=head2 GetBufferedPosNext
    About   : Buffer all vcf records by position
    Usage   : GetBufferedPosNext(vcf=$vcf;buffer=$buffer); assumes input of self like a hash and checks it the fields are present in this hash
    Args    :   my $targetvcf = Vcf->new(file=> $self -> {'targetvcf'})
	        $targetvcf->parse_header();
	        print $targetvcf->format_header();
	
		#iterate target file and collect data for one position across all files
	        my $posbuffer;
	
	        while ($posbuffer=GetBufferedPosNext('vcf' => $targetvcf,'buffer' => $posbuffer)){
			for my $x in @{$posbuffer-> {'current'}} {
				#do something
			}
		}

=cut

=head2 GetBufferedPosByLoc 
    About   : Get the next parsed vcf records in the given location. 
    Usage   : GetPositionBufferedPos(vcf=>$vcf,buffer=>$buffer,loc=[$chr,$pos]); assumes input of self like a hash and checks it the fields are present in this hash
    Args    : my $targetvcf = Vcf->new(file=> $self -> {'targetvcf'})
	        $targetvcf->parse_header();
	        print $targetvcf->format_header();
	
		#iterate target file and collect data for one position across all files
	        my $posbuffer;
		my $loc = {'CHROM' => 'Y','POS'=> 1000};
	        while ($posbuffer=GetBufferedPosByLoc('vcf' => $targetvcf,'buffer' => $posbuffer,$loc)){
			for my $x in @{$posbuffer-> {'current'}} {
				#do something
			}
		}

=cut

=head2 
    About   : LinePosLastHandleLogger
    Usage   :  LinePosLastHandleLogger();
    Args    :  LinePosLastHandleLogger(lineinterval => 1000,timeinterval=20). 
		Just uses time, localtime and $. to get what it needs. Checks every modulo lineinterval if time matches timeinterval and gives a warning.
=cut
