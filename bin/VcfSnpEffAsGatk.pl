#!/usr/bin/perl
#BEGIN{
#	push @INC,'/home/terpstramm/software/vcftools_0.1.9/lib/perl5/site_perl/';
#}
use warnings;
use strict;
use Vcf;
use Data::Dumper;
# input test
#$ARGV[0]='/home/terpstramm/Downloads/s12-193.annotated.vcf';
# use
my $use = <<"END";
	$0 in.vcf out.vcf
	The Gatk way report only worst annotation based on snpeff field in the input vcf and output with gatk like fields.
	Now with updated parser because GATKVariantAnnotator needs to be updated. For example Splice site annotations will fail to annotate.
END
die "no valid 'in.vcf' specified on command line\n$use" if(not(defined($ARGV[0])) ||not -e $ARGV[0]);
#main
AddZScoresBasedOnAdVals($ARGV[0],$ARGV[1]);

#old
##INFO=<ID=SNPEFF_AMINO_ACID_CHANGE,Number=1,Type=String,Description="Old/New amino acid for the highest-impact effect resulting from the current variant (in HGVS style)">
##INFO=<ID=SNPEFF_CODON_CHANGE,Number=1,Type=String,Description="Old/New codon for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_EFFECT,Number=1,Type=String,Description="The highest-impact effect resulting from the current variant (or one of the highest-impact effects, if there is a tie)">
##INFO=<ID=SNPEFF_EXON_ID,Number=1,Type=String,Description="Exon ID for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_FUNCTIONAL_CLASS,Number=1,Type=String,Description="Functional class of the highest-impact effect resulting from the current variant: [NONE, SILENT, MISSENSE, NONSENSE]">
##INFO=<ID=SNPEFF_GENE_BIOTYPE,Number=1,Type=String,Description="Gene biotype for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_GENE_NAME,Number=1,Type=String,Description="Gene name for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_IMPACT,Number=1,Type=String,Description="Impact of the highest-impact effect resulting from the current variant [MODIFIER, LOW, MODERATE, HIGH]">
##INFO=<ID=SNPEFF_TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID for the highest-impact effect resulting from the current variant">

#NEW
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">

sub AddZScoresBasedOnAdVals{
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	if($_[1] && not(-e $_[1])){
		my $vcfout = $_[1];
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn 'warning: no output vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	}
	
	
	my $vcf = Vcf->new(file=>$vcfin);
	$vcf->parse_header();
	my $annotations = GetAnnFunctionalAnnotations($vcf->get_header_line(key=>'INFO', ID=>'ANN'));
	SetHeadersSNPEFFANN('vcf' => $vcf, 'ann' => $annotations);
	print {$out} $vcf->format_header();
	#die "header check";
	#$vcf->recalc_ac_an(0);
	while (my $x=$vcf->next_data_hash()){
		
		#die Dumper($x);
		#my $ann = GetAnn($x);
		#die $ann;
		FillSNPEFFANNFields('record' => $x,'ann' => $annotations);
		
		#warn 
		#if(not($x->{'INFO'}{'VariantType'} =~ m/^MULTIALLELIC/)){
		#	$vcf->add_format_field($x,'Z'); 
		#	$vcf->add_format_field($x,'Y'); 
		#	my @samples = keys(%{$x->{gtypes}});
		#	#warn Dumper(@samples); 
		#	for my $sample (@samples){
		#		my @AD = split(',',$x->{'gtypes'}{$sample}{'AD'})if(defined($x->{'gtypes'}{$sample}{'AD'})&& not($x->{'gtypes'}{$sample}{'AD'} =~ m/^\.,\./));
		#		my $Zscore;
		#		if( defined($AD[0]) && defined($AD[1]) && scalar(@AD)== 2 && $AD[0] + $AD[1] > 0 ){
		#			$Zscore = 0.5 * ($AD[1] / ($AD[1] + $AD[0]) - 0.5) * sqrt($AD[1] + $AD[0]);
		#		}else{
		#			$Zscore = 0;
		#		}
		#		$x->{'gtypes'}{$sample}{'Z'}=$Zscore;
		#		$x->{'gtypes'}{$sample}{'Y'}=(1-$Zscore);
		#		
		#		my $F;
		#		if( defined($AD[0]) && defined($AD[1]) && scalar(@AD)== 2 && $AD[0] + $AD[1] > 0 && $AD[1] > 0 ){
		#			
		#			$F = ($AD[1] / sum(@AD));
		#			
		#		}else{
		#			
		#			$F = '0';
		#		}
		#		
		#		$x->{'gtypes'}{$sample}{'F'}=$F;
		#	}
		#}
		#}else{
		#	warn Dumper($x);die;
		#}
		
		#warn Dumper($vcf);
		#warn Dumper($x);
		
		print {$out} $vcf->format_line($x);
		#die;
		#die if ($. > 50);
	}
	$vcf->close();
}

sub GetAnnFunctionalAnnotations {
	my $annParsed = shift @_ or die "no input given for subroutine";
	die Dumper($annParsed) if(not(defined($annParsed -> [0] -> {'Description'})));
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
	#die Dumper(\$anndescr,$annotations). $anndescr;
	#should contain:
	#$VAR2 = [
    #      'ALLELE',
    #      'ANNOTATION',
    #      'ANNOTATION_IMPACT',
    #      'GENE_NAME',
    #      'GENE_ID',
    #      'FEATURE_TYPE',
    #      'FEATURE_ID',
    #      'TRANSCRIPT_BIOTYPE',
    #      'RANK',
    #      'HGVS_C',
    #      'HGVS_P',
    #      'CDNA_POS_AND_CDNA_LENGTH',
    #      'CDS_POS_AND_CDS_LENGTH',
    #      'AA_POS_AND_AA_LENGTH',
    #      'DISTANCE',
    #      'ERRORS_AND_WARNINGS_AND_INFO'
    #    ];
    return $annotations
}
sub SetHeadersSNPEFFANN {
	my $self;
	%{$self}= @_;
	my $vcf = $self -> {'vcf'} or die "no vcf object as input";
	my $annotations = $self -> {'ann'} or die "no ann object(annotations) as input";
	for my $annotation (@{$annotations}){
		
		$vcf -> add_header_line({'key' => 'INFO', 'ID' => SnpEffAnnField($annotation -> {'name'}), 'Number' => 'A','Type' => 'String','Description' => 'Snpeff field: ' . SnpEffAnnField($annotation -> {'alias'}) . '. This is an subfield of the SnpEff annotation selected by ' . $0 . ' on the most harmful prediction. The format is described on http://snpeff.sourceforge.net/SnpEff_manual.html .'});
		
	}
	#die Dumper($vcf);
	#return $vcf;
	
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
	print Dumper($worstPredsByAlleleParsed,$r);
	die Dumper($worstPredsByAlleleParsed,$r);
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
	die Dumper(@worstAnnotations)."\n".$. if(scalar(@{$alts})>1);
	
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