#!/usr/bin/env perl -w
use warnings;
use strict;
use Vcf;
#use Devel::Confess;
use Data::Dumper;

my $use = <<"END";
	$0 COMPLEXMERGEVCF LOSTSAMPLEINFOVCF [VCF]
		Tries to recover sample annotations (e.g. GT, GQ, PL, AD, etc) lost by using GATK CombineVariants from one or more VCF files prioritised by input order.

	Needs vcftools for library VCF.pm and htslib for tabix and bgzip.
	Gives prio to the first VCF present in the list for annotating
END

my $complexMergeVcf = shift @ARGV or die "No COMPLEXMERGEVCF option supplied. ARGV=[@ARGV]\n$use";
die localtime(time())." [ERROR] $0: COMPLEXMERGEVCF already exists pls remove '$complexMergeVcf'." if( -e $complexMergeVcf);  
my $lostSampleVcf = shift @ARGV or die localtime(time())." [ERROR] $0: No LOSTSAMPLEINFOVCF option supplied. ARGV=[@ARGV]\n$use";

die "no valid 'VCF' specified on command line\n$use VCF=[@ARGV]" if( scalar(@ARGV)==0 || not( -e $ARGV[0]));
#main
main($complexMergeVcf, $lostSampleVcf, @ARGV);

sub main{
	my $vcfcomplex = shift @_;
	my $vcfin =shift @_ or die "$use";#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	die "LOSTSAMPLEINFOVCF does not exist. Arg:'$vcfin' " if(not(-e $vcfin));
	my $out =*STDOUT;
	
	warn localtime(time())." [INFO] $0: Parsing header of 'vcfin'";
	
	my $resource;
	die localtime(time())." [ERROR] $0: No resource vcf given.\n$use" if (scalar(@ARGV)== 0 );
	@{$resource -> {'invcfs'}} = @_;
	my @idxVcfs=();
	@{$resource -> {'idxinvcfs'}}=();
	for my $rvcf (@{$resource -> {'invcfs'}}) {
		my $tmpIdx = scalar(@{$resource -> {'idxinvcfs'}});
		my $idxinVcf = $vcfin.'tmp'.$tmpIdx.'.vcf.gz';
		CreateTabixIndexedVcf('invcf' => $rvcf,'outvcf' => $idxinVcf);
		push(@{$resource -> {'idxinvcfs'}}, $idxinVcf);
		my $idxVcf = Vcf->new(file=>$idxinVcf);
		$idxVcf->parse_header();
		push(@{$resource -> {'idxvcfs'}}, $idxVcf);
	}

	my $vcf = Vcf->new(file=>$vcfin);
	$vcf->parse_header();
	print {$out} $vcf->format_header();
	
	open(my $vcfcomplexout ,'>', $vcfcomplex) or die die localtime(time())." [ERROR] $0: Cannot open vcfcomplex:'$vcfcomplex'.";
		print {$vcfcomplexout} $vcf->format_header();

	warn localtime(time())." [INFO] $0: Iterating records."; my $records=0;
	while (my $x=$vcf->next_data_hash()){
		#die Dumper $x;
		#annotate populate fill the
		#warn "ret=".SampleInfoUnannotated($x).". line=$."; 
		if(SampleInfoUnannotated($x)){
			my @vcfline = split("\t",$vcf->format_line($x));
			 $vcfline[6] = "PASS";
			print {$vcfcomplexout} join("\t", @vcfline) if(AnnotateSampleInfo($x, $resource) == 0);
		}
		print {$out} $vcf->format_line($x);
		#die Dumper $x;
		$records++;
	}
	$vcf->close();
	warn localtime(time())." [INFO] $0: Done. Processed $records";
}

sub CalleriseInfoHeader {
	my $self;
        %{$self}= @_;
        my $vcf = $self -> {'vcf'} or die "no vcf object as input";
	my $caller = $self -> {'caller'} or die "no caller string as input";
	for my $field (keys(%{$vcf -> {header} -> {INFO}})){
		$vcf -> add_header_line({'key' => 'INFO',
			'ID' => $caller.$field,
			'Number' => $vcf -> {header} -> {INFO} -> {$field} -> {'Number'},
			'Type' => $vcf -> {header} -> {INFO} -> {$field} -> {'Type'} ,
			'Description' => $vcf -> {header} -> {INFO} -> {$field} -> {'Description'}});
	}
}
sub CalleriseFormatHeader {
        my $self;
        %{$self}= @_;
        my $vcf = $self -> {'vcf'} or die "no vcf object as input";
        my $caller = $self -> {'caller'} or die "no caller string as input";
       	for my $field (keys(%{$vcf -> {header} -> {FORMAT}})){
                $vcf -> add_header_line({'key' => 'FORMAT',
                        'ID' => $caller.$field,
                        'Number' => $vcf -> {header} -> {FORMAT} -> {$field} -> {'Number'},
                        'Type' => $vcf -> {header} -> {FORMAT} -> {$field} -> {'Type'} ,
                       	'Description' => $vcf -> {header} -> {FORMAT} -> {$field} -> {'Description'}});
        }
}
sub CalleriseRecordInfoData {
	my $self;
	%{$self} = @_;
	my $record = $self -> {'record'} or die "no record object as input";
	my $caller = $self -> {'caller'} or die "no caller string as input";
	
	for my $field (keys(%{$record -> {'INFO'}})){
		$record -> {'INFO'} -> {$caller.$field} = $record -> {'INFO'} -> {$field};
	}
}
sub CalleriseRecordFormatData {
        my $self;
        %{$self} = @_;
        my $record = $self -> {'record'} or die "no record object as input";
        my $caller = $self -> {'caller'} or die "no caller string as input";
	
	my @callerisedFormatFields;
        for my $field (@{$record -> {'FORMAT'}}){
		push(@callerisedFormatFields,$caller.$field);
        }
	push(@{$record -> {'FORMAT'}},@callerisedFormatFields);
	
	for my $sample (keys(%{$record -> {'gtypes'}})){
		for my $field (keys(%{$record -> {'gtypes'} -> {$sample}})){
			$record -> {'gtypes'} -> {$sample} -> {$field} = './.' if($field eq "GT" && $record -> {'gtypes'} -> {$sample} -> {$field} eq '.');
			next if($field eq 'GT' && ($record -> {'gtypes'} -> {$sample} -> {$field} eq './.'));
			$record -> {'gtypes'} -> {$sample} -> {$caller.$field} = $record -> {'gtypes'} -> {$sample} -> {$field};
		}

	}
}


sub GetAnnFunctionalAnnotations {
	my $annParsed = shift @_ or die "no input given for subroutine";
	die "invalid input".Dumper($annParsed) if(not(defined($annParsed -> [0] -> {'Description'})));
	my $anndescr = $annParsed -> [0] -> {'Description'};
	$anndescr =~ s/Functional\ annotations\:\ \'|\'//g;
	my $anndescrvcf = $anndescr;
	$anndescrvcf =~ s/ \/ /_AND_/g;
	$anndescrvcf =~ s/\./_/g;
	
	my $annotations;
	#@{$annotations}
	
	my @anns = split(' \| ', $anndescr);
	my @annsVcf = split(' \| ', $anndescrvcf);
	for my $annVcf (@annsVcf){
		my $ann = shift(@anns);
		my %h = ('name' => uc($annVcf),'alias' => $ann);
		push(@{$annotations},\%h);
	}

    return $annotations
}
sub SetHeadersSNPEFFANN {
	my $self;
	%{$self}= @_;
	my $vcf = $self -> {'vcf'} or die "no vcf object as input";
	my $annotations = $self -> {'ann'} or die "no ann object(annotations) as input";
	for my $annotation (@{$annotations}){
		
		$vcf -> add_header_line({'key' => 'INFO', 
			'ID' => SnpEffAnnField($annotation -> {'name'}),
			'Number' => 'A','Type' => 'String',
			'Description' => 'Snpeff field: ' . SnpEffAnnField($annotation -> {'alias'}) . '. This is an subfield of the SnpEff annotation selected by ' . $0 . ' on the most harmful prediction. The format is described on http://snpeff.sourceforge.net/SnpEff_manual.html .'});
		
	}
	
}
sub SnpEffAnnField{
	return 'SNPEFFANN_'.shift @_;
}
sub GetAnn{
	my $r = shift @_;
	if(defined($r->{'INFO'}{'ANN'})){
		return $r->{'INFO'}{'ANN'};
	}
	die "Vcf record doesn't have ANN INFO Field or is strangely annotated!".Dumper($r);
}
sub GetAlt{
	my $r = shift @_;
	if(defined($r->{'ALT'})){
		return $r->{'ALT'};
	}
	die "Vcf record doesn't have ALT Field or is strangely annotated!".Dumper($r);
}
sub FillSNPEFFANNFields {
	my $self;
	%{$self}= @_;
	my $r = $self -> {'record'} or die "no record object as input";
	my $annotations = $self -> {'ann'} or die "no ann object(annotations) as input";
	AnnotateRecordAnnFields('record' => $r,'ann' => $annotations);
}
sub AnnotateRecordAnnFields {
	my $self;
	%{$self}= @_;
	my $r = $self -> {'record'} or die "no record object as input";
	my $annotations = $self -> {'ann'} or die "no ann object(annotations) as input";
	my $worstPredsByAllele = GetWorstPredictionsByAllele('record' => $r,'ann' => $annotations);
	my $worstPredsByAlleleParsed = ParseAnnotation('annfields'=>$worstPredsByAllele,'ann'=>$annotations);
	
	for my $field (@{$worstPredsByAlleleParsed}){
		$r -> {'INFO'} -> {SnpEffAnnField($field -> {'name'})}=join(',',@{$field -> {'values'}});
	}
	#print Dumper($worstPredsByAlleleParsed,$r);
	#die Dumper($worstPredsByAlleleParsed,$r);
}
sub GetWorstPredictionsByAllele {
	my $self;
	%{$self}= @_;
	my $r = $self -> {'record'} or die "no record object as input";
	my $annotations = $self -> {'ann'} or die "no ann object(annotations) as input";
	
	#$worstPredsByAlt -> {alleles} = [ALT]
	#$worstPredsByAlt -> {ANNFIELDS} -> [FORALT results];
	my $alts = GetAlt($r);
	my @worstAnnotations;
	for my $alt (@{$alts}){
		my $worstAnnotation = GetWorstAnnByAllele('record'=>$r,'alt'=>$alt,'ann'=>$annotations);
		push(@worstAnnotations,$worstAnnotation );
		die "Missing annotation field in record?!" if(not(defined($worstAnnotation)));
	}
	
	return \@worstAnnotations;
}

sub GetWorstAnnByAllele {
	my $self;
	%{$self}= @_;
	my $r = $self -> {'record'} or die "no record object as input";
	my $alt = $self -> {'alt'} or die "no alt value as input";
	my $annotations = $self -> {'ann'} or die "no annotations object value as input";
	
	my @anns = split(",",GetAnn($r));
	
	my $worstAnn;
	
	while(not(defined($worstAnn))){
		my $ann = shift @anns;
		$worstAnn = $ann if (GetAnnotationAllele('annfields' => [$ann],'ann' => $annotations) eq $alt);
	}
	return $worstAnn;
	
}

sub ParseAnnotation {
	my $self;
	%{$self}= @_;
	my $annFields = $self -> {'annfields'} or die "no annotations object value as input";
	my $annotationHeader = $self -> {'ann'} or die "no annotations object value as input";

	my $parsed;

	for my $annField (@{$annFields}){
		my @anns = split /\|/,$annField,-1;
		die "Annotations not equal to header:\n".Dumper(\$annField,\@anns,$annotationHeader) if(not(scalar(@anns)==scalar(@{$annotationHeader})));
		
		my $fieldindex = 0;
		#die Dumper($annotationHeader);
		while($fieldindex < scalar(@{$annotationHeader})){
			$parsed -> [$fieldindex] -> {'name'} = $annotationHeader -> [$fieldindex] -> {'name'};
			push(@{$parsed -> [$fieldindex] -> {'values'}},$anns[$fieldindex]);
			$fieldindex++;
		}
	}
	return $parsed;
}
sub GetAnnotationAllele {
	my $self;
	%{$self}= @_;
	my $annField = @{$self -> {'annfields'}}[0] or die "no annotations object value as input";
	my $annotationHeader = $self -> {'ann'} or die "no annotations object value as input";
	my $parsed = ParseAnnotation(%{$self});
	for my $annotation (@{$parsed}){
		return $annotation -> {'values'} -> [0] if($annotation -> {'name'} eq 'ALLELE' && scalar(@{$annotation -> {'values'}}) == 1);
	}
}

sub CreateTabixIndexedVcf {
	my $self;
	%{$self}= @_;
	my $invcf = $self -> {'invcf'} or die "no vcf file as input";
	my $outvcf = $self -> {'outvcf'} or die "no bgzipped file as output";
	die "Wrong extension (not .vcf.gz) for outvcf $outvcf " if(not($outvcf =~ m/\.vcf\.gz$/));
	
	my $indexcmd = "set -x && bgzip -c $invcf > $outvcf && tabix -p vcf $outvcf";
	my $ret = `$indexcmd`;
        $? == 0 or die "system cmd `$indexcmd` failed: $ret\n exit code: $?";
	

}
sub SampleInfoUnannotated {
	my $record = shift @_ or die "No input record";
	
	if((scalar(@{$record -> {'FORMAT'}}) == 1 && @{$record -> {'FORMAT'}}[0] eq 'GT')|| 1){
		my $hasGenotypes = 0;
		my $hasAD = 0;
		for my $sample (keys(%{$record -> {'gtypes'}})){
                	for my $field (keys(%{$record -> {'gtypes'} -> {$sample}})){
				if($field eq 'GT' && 
					($record -> {'gtypes'} -> {$sample} -> {$field} eq '.' || 
					 $record -> {'gtypes'} -> {$sample} -> {$field} eq './.')){
					$hasGenotypes ++;
				}
				$hasAD++ if($field eq 'AD' &&
					not(defined($record -> {'gtypes'} -> {$sample} -> {$field})) ); 
			}
		}
		#warn "$. GT $hasGenotypes AD $hasAD"; 
		if( $hasAD == 0){
			#no Genotypes or AD present so return true
			return 1;
		}
	}
	#all genotypes present so return false
	return 0;
}

sub AnnotateSampleInfo {
	my $record = shift @_;
	my $resource = shift @_;
	
	#@{$resource -> {'idxvcfs'}}
	#$idxvcf->open('region' => $record -> {'CHROM'}.':'.$record -> {'POS'}.'-'.$record -> {'POS'} ); #$vcf->open(region=>'1:12345-92345')
	my $annotated = 0;
	# go through all resourcevcfs && annotate if found at same pos / ref/ alts

	for my $idxvcf  (@{$resource -> {'idxvcfs'}}){
		$idxvcf->open('region' => $record -> {'CHROM'}.':'.$record -> {'POS'}.'-'.$record -> {'POS'} );
		while (my $x=$idxvcf->next_data_hash()){
			#die Dumper($x,$record)." ";
			if($record -> {'POS'} eq $x  -> {'POS'} && $record -> {'REF'} eq $x  -> {'REF'} && scalar(@{$record -> {'ALT'}}) == scalar(@{$x -> {'ALT'}})){
				my $alteq = 1;
				my $altidx = 0;
				while($alteq == 1 && $altidx < scalar(@{$record -> {'ALT'}}) ){
					#Check if alts are equal
					$alteq = 0 if $record -> {'ALT'} -> [$altidx] ne $x -> {'ALT'} -> [$altidx];
					#if($annotated == 1){
					#	return 1;
					#}
					$altidx++;
				}
				#if chrom, pos, ref,alts are equal then for each field in format liftover $x to $record if not already present... Hope this works
				if($alteq == 1){
					$annotated = 1;
					for my $sample (keys(%{$record -> {'gtypes'}})){
						for my $field (keys(%{$x -> {'gtypes'} -> {$sample}})){
							if($field eq 'GT' && 
								($record -> {'gtypes'} -> {$sample} -> {$field} eq '.' ||
								 $record -> {'gtypes'} -> {$sample} -> {$field} eq './.') ){
								$record -> {'gtypes'} -> {$sample} -> {$field} = $x -> {'gtypes'} -> {$sample} -> {$field};
							}
							if( not( defined($record -> {'gtypes'} -> {$sample} -> {$field})) ){
								$record -> {'gtypes'} -> {$sample} -> {$field} = $x -> {'gtypes'} -> {$sample} -> {$field};
							}
			                	}
					}
					#die "Should work".Dumper($record);
					#
				}
			}
		}
	}
	if(not($annotated == 1)){
		return 0;
	}
	return 1;
}
