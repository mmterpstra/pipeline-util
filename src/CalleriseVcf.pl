#!/usr/bin/env perl -w
use warnings;
use strict;
use Vcf;
use Data::Dumper;

my $use = <<"END";
	$0 caller in.vcf out.vcf
		add more annotations for vcfs by prepending existing records with the caller text 
END

my $caller = shift @ARGV or die " No caller option supplied. $use";

die "no valid 'in.vcf' specified on command line\n$use" if(not(defined($ARGV[0])) ||not -e $ARGV[0]);
#main
main($ARGV[0],$ARGV[1],$caller);

sub main{
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	if($_[1] && not(-e $_[1])){
		my $vcfout = $_[1];
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn localtime(time())." [WARN] $0: no output vcf file specified printing to STDOUT"."\n";
		$out =*STDOUT;
	}
	my $caller = $_[2];
	
	warn localtime(time())." [INFO] $0: Parsing header.";
	my $vcf = Vcf->new(file=>$vcfin);
	#warn "parse header";
	$vcf->parse_header();
	
	#warn "get annotations";
	CalleriseInfoHeader('vcf' => $vcf, 'caller' => $caller);
	CalleriseFormatHeader('vcf' => $vcf, 'caller' => $caller);
	print {$out} $vcf->format_header();

	#die "header check";
	#$vcf->recalc_ac_an(0);
	warn localtime(time())." [INFO] $0: Iterating records."; my $records=0;
	while (my $x=$vcf->next_data_hash()){
#		if(defined($annheaderline) && scalar(@{$annheaderline})){
#			FillSNPEFFANNFields('record' => $x,'ann' => $annotations);
#		}
		CalleriseRecordInfoData('record' => $x,'caller' => $caller);
		CalleriseRecordFormatData('record' => $x,'caller' => $caller);
		#die Dumper $x;
		print {$out} $vcf->format_line($x);
		#die Dumper $x;
		$records++;
#               #if(scalar(@{GetAlt($x)})>1){
#		#		die Dumper($x)."\n".GetAlt($x)."\n".$. ;
#		#}

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
