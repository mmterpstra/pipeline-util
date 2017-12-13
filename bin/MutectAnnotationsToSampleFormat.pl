#!/usr/bin/perl
use warnings;
use strict;
use Vcf;
use Data::Dumper;

my $use = <<"END";
	$0 fields in.vcf out.vcf
		move annotations from info field to genotype field. Fields should be comma seperated. 
		Also tries to move the filter to genotype
		minimal test / recommended usage:
			perl ../pipeline-util/bin/MutectAnnotationsToSampleFormat.pl TLOD,NLOD,MIN_ED,MAX_ED,ECNT,HCNT,PON MuTect2.out.vcf 
END

my $fields = shift @ARGV or die " No fields option supplied. $use";

die "no valid 'in.vcf' specified on command line\n$use" if(not(defined($ARGV[0])) ||not -e $ARGV[0]);
#main
main($ARGV[0],$ARGV[1],$fields);

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
	my $fieldsCs = $_[2];
	
	
	warn localtime(time())." [INFO] $0: Parsing header.";
	my $vcf = Vcf->new(file=>$vcfin);
	#warn "parse header";
	$vcf->parse_header();
	
	#mutect should only contain single sample so get genotype or die with error
	my $gt  = GetGenotypeOrDie($vcf);

	#interpret fields from cmdline
	my $fields;
	@{$fields} = split(",",$fieldsCs);
	#@{$fields}=("TLOD","NLOD");
	
	#Update headers for FORMAT info and also add FT annotation and print header as vcf
	InfoFieldsToGenotypeFieldsInHeader('vcf' => $vcf, 'fields' => $fields);
	print {$out} $vcf->format_header();
	#die "header check";
	
	#iterate file
	warn localtime(time())." [INFO] $0: Iterating records."; my $records=0;
	while (my $x=$vcf->next_data_hash()){
		InfoFieldsToGenotypeFieldsInRecord('record' => $x,'fields' => $fields, 'gt' => $gt);
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

sub InfoFieldsToGenotypeFieldsInHeader {
	my $self;
        %{$self}= @_;

	#if self FT?
	$self->{'vcf'}->add_header_line({
            'default' => '.',
            'ID' => 'FT',
            'Description' => "Filter based on mutect filters: alt_allele_in_normal='Evidence seen in the normal sample'|clustered_events='Clustered events observed in the tumor'|clustered_read_position='Evidence for somatic variant clusters near the ends of reads'|germline_risk='Evidence indicates this site is germline, not somatic'|homologous_mapping_event='More than three events were observed in the tumor'|multi_event_alt_allele_in_normal='Multiple events observed in tumor and normal'|panel_of_normals='Seen in at least 2 samples in the panel of normals'|str_contraction='Site filtered due to contraction of short tandem repeat region'|strand_artifact='Evidence for alt allele comes from one read direction only'|t_lod_fstar='Tumor does not meet likelihood threshold'|triallelic_site='Site filtered because more than two alt alleles pass tumor LOD'|",
            'Number' => '1',
            'handler' => undef,
            'key' => 'FORMAT',
            'Type' => 'String'});
	SelfRequire(%{$self}, 'req'=> ['vcf', 'fields']);
	
	for my $field ( @{ $self -> {'fields'} }){
		my $h=$self->{'vcf'} -> get_header_line(key=>'INFO', ID=>$field);
		$h -> [0] -> {'key'} ='FORMAT';
		$self->{'vcf'}->add_header_line($h -> [0]);
		$self->{'vcf'}->remove_header_line(key=>'INFO', ID=>$field);
		#warn Dumper($h,$self->{'vcf'}->format_header())." ";
		
	}
	#$vcf->get_header_line(key=>'INFO', ID=>'AC');
	#$vcf->add_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'});
	#$vcf->remove_header_line(key=>'INFO', ID=>'AC');
}
sub InfoFieldsToGenotypeFieldsInRecord {
	my $self;
        %{$self}= @_;
		
	SelfRequire(%{$self}, 'req'=> ['record','fields','gt']);

	$self -> {'record'} -> {'gtypes'} -> { $self -> {'gt'} } -> {'FT'} = join(';',@{$self -> {'record'} -> {'FILTER'}});
	
	for my $field ( @{ $self -> {'fields'} }){
		#mv $self -> {'record'} -> {'INFO'} -> {$field}
		#to $self -> {'record'} -> {'gtypes'} -> { $self -> {'gt'} } -> {$field}
		if($self -> {'record'} -> {'INFO'} -> {$field}){
			$self -> {'record'} -> {'gtypes'} -> { $self -> {'gt'} } -> {$field} = $self -> {'record'} -> {'INFO'} -> {$field};
			delete($self -> {'record'} -> {'INFO'} -> {$field});
		}
	}
	#die Dumper($self -> {'record'})." ";
}
sub SelfRequire {
	my $self;
        %{$self}= @_;
	#die Dumper($self)." ";
	my $reqs = $self -> {'req'} or die "no requirement array as input";
	for my $tag (@{$reqs}){
		defined($self -> {$tag}) or die "no requirement '$tag' as input";
	}
}

sub GetGenotypeOrDie {
	#my $self;
        #%{$self}= @_;
	my $vcf = shift(@_) or die "no vcf object";
	my (@samples) = $vcf->get_samples();
	die "This is not a mutect2 vcf because it contains no samples" if(scalar(@samples) == 0);
	die "This is not a mutect2 vcf because it contains more than 1 sample" if(scalar(@samples) > 1);
	my $gt =$samples[0];
	#die Dumper(\@samples)." ";
	return $gt;
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
