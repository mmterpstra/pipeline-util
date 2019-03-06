package pipeline::util;

use 5.020002;
use strict;
use warnings;
use Vcf;
use Data::Dumper;
use Carp qw(verbose confess carp croak);
require Exporter;

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
	FormatWalkTargetLineAsVcfLine
	FormatWalkTargetLineAsVcfHeader
	_formatwalkasvcflineswithfile
);

our $VERSION = "0.8.9";


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
			warn '60606060';
			$vcf -> {buffer} = GetBufferedPosByLoc('vcf' => \$vcf -> {handle},'buffer' => $vcf -> {buffer}, 'loc'=>  $targetposbuffer -> {'current'} -> [0]);
			#do some annotations
			warn 'target='.
				$targetposbuffer -> {'current'} -> [0] -> {'CHROM'}.':'.
				$targetposbuffer -> {'current'} -> [0] -> {'POS'}.
				" vcf=".$vcf -> {'buffer'} -> {'current'} -> [0] -> {'CHROM'}.':'.
				$vcf -> {'buffer'} -> {'current'} -> [0] -> {'POS'}
		}
		 
	}
	#dump rest of buffer
	#insert code

}
sub NewWalk{
	#run like 	$walkdata = NewWalk('targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}});
	# this opens the relavant files and stores the data 
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['targetvcf', 'vcflist']);
	my $walk;
        $walk -> {'targetvcf' } -> {'file' } = $self -> {'targetvcf'};
	$walk -> {'targetvcf' } -> {'handle' } = Vcf->new(file=> $self -> {'targetvcf'});
	$walk -> {'targetvcf' } -> {'handle' } -> parse_header();
	$walk -> {'targetvcf' } -> {'buffer' } = {'current' => [],'next' => [ $walk -> {'targetvcf' } -> {handle} -> next_data_hash() ] };#GetBufferedPosNext('vcf' =>\$targetvcf,'buffer' => $targetposbuffer) or die "Empty vcf!! Cannot parse";

	#$walk -> {'vcfs' };
	for my $vcf (@{$self -> {'vcflist'}}){
		push(@{$walk -> {'vcfs' }},{'file' => $vcf,'handle'=>Vcf->new(file=> $vcf)});
	        $walk -> {'vcfs' } -> [-1] -> {handle}->parse_header();
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
			#	$vcf -> {'handle' } -> format_line($record);
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
			for my $infofield (keys(%{$self -> {'record'} ->  {'INFO'}})){
				if(not(defined $self -> {'targetrecord'} -> {'INFO'} -> {$infofield}) || $self -> {'targetrecord'} -> {'INFO'} -> {$infofield} eq '.'){
					$self -> {'targetrecord'} -> {'INFO'} -> {$infofield} = $self -> {'record'} -> {'INFO'} -> {$infofield};
				}
				
			}			
			for my $formatfield (keys(%{$self -> {'record'} ->  {'gtypes'} -> {$sample} })){
				#warn Dumper(keys(%{$self -> {'record'} ->  {'gtypes'} -> {$sample} }))." ";

				if(defined($formatfield) && $formatfield eq 'GT' && 
					($self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} eq '.' ||
					 $self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} eq './.') ){
					$self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} = $self -> {'record'} ->  {'gtypes'} -> {$sample} -> {$formatfield};
				}
				if( not( defined($self -> {'targetrecord'} ->  {'gtypes'} -> {$sample} -> {$formatfield})) || $self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} eq '.' ){
					$self -> {'targetvcfhandle'} -> add_format_field($self -> {'targetrecord'},$formatfield);
					$self -> {'targetrecord'} -> {'gtypes'} -> {$sample} -> {$formatfield} = $self -> {'record'} -> {'gtypes'} -> {$sample} -> {$formatfield};
					#$vcf->add_format_field($x,'FOO');
				}
			}
		}
	}
	return $self -> {'targetrecord'};
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
			#Dumper($record);
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
	
	$out .= $walk -> {'targetvcf' } -> {'handle' } -> format_header();
	
	
	return $out;
	
}

sub SyncVcfToTarget{
	my $self ;%{$self}= @_ ;
	#@{$self -> {'vcflist'}}=($self -> {'targetvcf'});
#	SelfRequire(%{$self},'req'=> ['walk']);
	SelfRequire(%{$self},'req'=> ['vcf','pos']);	
	#my $walk = $self -> {'walk'};
	my $continue = 1;
	while ($continue == 1 && ($self -> {'vcf' } -> {'buffer'} = GetBufferedPosNext('vcf' =>$self -> {'vcf' } -> {'handle'},'buffer' => $self -> {'vcf' } -> {'buffer'}))){
		
		#Goes until equal position or greater then position; #not(ChromPosIsNext('loc1' => {'CHROM'=>4,'POS'=>55593536}, 'loc2' => $targetposbuffer -> {'current'} -> [0], 'vcf' => $targetvcf))
		
		#Goes till
		if(not(defined($self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0])) || ChromPosIsNext('loc2' => {$self -> {'CHROM'},$self -> {'POS'}}, 'loc1' => $self -> {'vcf'} -> {'buffer'} -> {'next'} -> [0], 'vcf' => $self -> {'vcf'} -> {'handle'})){
			
			warn Dumper({'CHROM' => $self -> {'vcf'} -> {'buffer'} -> {'current'} -> [0] -> {'CHROM'},'POS' => $self -> {'vcf'} -> {'buffer'} -> {'current'} -> [0] -> {'POS'}})." ";
			
			$continue = 0;
			
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
		if(ChromPosIsNext('loc2' => $self -> {'pos'}, 'loc1' => $targetposbuffer -> {'next'} -> [0], 'vcf' => $targetvcf)){
			
			warn Dumper({'CHROM' => $targetposbuffer -> {'current'} -> [0] -> {'CHROM'},'POS' => $targetposbuffer -> {'current'} -> [0] -> {'POS'}}) . "  ";
			
			$continue = 0;
			
		}
		
	}
	warn Dumper(GetLoc(%{$targetposbuffer -> {'current'} -> [0]}));
	warn "done ";
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
sub GetLoc {
	my $self;
        %{$self}=@_ or  confess "input not parseable as hash:".Dumper(@_). ' ';

	return {'CHROM' => $self -> {'CHROM'},'POS' => $self -> {'POS'}};
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
		#warn '## VC ## '. Dumper($x). ' ';
		$self -> {'buffer'} -> {'next'} = [$x];
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
		push(@{$self -> {'buffer'} -> {'current'}},%{$self -> {'buffer'} -> {'next'}});
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
	
	my $contigcur=0;
	my $contigidx1=-1;my $contigidx2=-0.5;
	#warn Dumper($self -> { 'vcf'}). ' ' if($. > 300);
	for my $contig (GetContigsAsArray('vcf' => $self -> { 'vcf'})){
		$contigidx1=$contigcur if(defined $self -> { 'loc1'} -> {'CHROM'} && $self -> { 'loc1'} -> {'CHROM'} eq $contig);
		$contigidx2=$contigcur if(defined $self -> { 'loc2'} -> {'CHROM'} && $self -> { 'loc2'} -> {'CHROM'} eq $contig);
		$contigcur++;
	}
	warn "############ $contigidx1 > $contigidx2";
	return 1 if($contigidx1 > $contigidx2 );

	return 0;
}

sub ChromPosIsEq {
	my $self;%{$self} = @_;
	warn Dumper({'loc1' =>  GetLoc($self -> {'loc1'}), 'loc2' => GetLoc($self -> {'loc2'})}). " ";
	if(defined($self -> { 'loc2'}) &&
		defined($self -> { 'loc2'} -> {'CHROM'}) &&
		$self -> { 'loc1'} -> {'CHROM'} eq $self -> { 'loc2'} -> {'CHROM'} && 
		$self -> { 'loc1'} -> {'POS'} eq $self -> { 'loc2'} -> {'POS'}){
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
	#carp "loc1".Dumper($self -> { 'loc1'})."loc2".Dumper($self -> { 'loc2'});
	if(defined( $self -> { 'loc2'})&& defined( $self -> { 'loc2'} -> {'CHROM'}) && ($self -> { 'loc1'} -> {'CHROM'} eq $self -> { 'loc2'} -> {'CHROM'} && 
			$self -> { 'loc1'} -> {'POS'} > $self -> { 'loc2'} -> {'POS'}) || 
		($self -> { 'loc1'} -> {'CHROM'} ne $self -> { 'loc2'} -> {'CHROM'} && ContigIsNext('loc1'=> $self -> { 'loc1'},'loc2' => $self -> { 'loc2'}, 'vcf' => $self -> {'vcf'})) ){
		warn 'ChromPosIsNext::ret='.1;
		return 1;
	}else{
		return 0;
	}
}

sub ChromPosIsNextOrEqual  {
	my $self; %{$self} = @_;
	#carp "loc1".Dumper($self -> { 'loc1'})."loc2".Dumper($self -> { 'loc2'});
	if(ChromPosIsNext(%{$self}) || ChromPosIsEq(%{$self})){
		warn 'ChromPosIsNext::ret='.1;
		return 1;
	}else{
		return 0;
	}
}
sub LocGetChromPosAsString{
	my $self; %{$self} = @_;
	return $self -> {'loc'} -> {'CHROM'}.":".$self -> {'loc'} -> {'POS'};
}

sub GetLoc {
	my $self; $self = shift @_;
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
