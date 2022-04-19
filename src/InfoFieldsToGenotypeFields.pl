#!/usr/bin/env perl
use warnings;
use strict;
use Vcf;
use Data::Dumper;

my $use = <<"END";
	$0 -f 'AD,SB,DP4=DP' -g 'SAMPLE' -i in.vcf -o out.vcf
		sends sample specific annotations to the genotype field 
END

#my $caller = shift @ARGV or die " No caller option supplied. $use";

#die "no valid 'in.vcf' specified on command line\n$use" if(not(defined($ARGV[0])) ||not -e $ARGV[0]);
#main
use  Getopt::Long;
my $opts ={};
my $result = GetOptions ("fields=s"	=> \$opts -> {'f'},
                        "input=s"	=> \$opts -> {'i'},
			"output=s"	=> \$opts -> {'o'},
			"genotype=s"	=> \$opts -> {'g'},
);
warn Dumper($result, $opts)." ";

#getopts('Ci:o:c:', $opts);

warn Dumper($opts)." ";
die "no valid vcf'".$opts -> {'i'}."' specified on command line\n$use" if(not(defined($opts -> {'i'})) || not -e $opts -> {'i'});

main($opts);

sub main{
	my $opts = $_[0];
	my $vcfin = $opts -> {'i'};#'some.annotated.vcf';
	my $out;
	if($opts -> {'o'} && not(-e $opts -> {'o'})){
		my $vcfout = $opts -> {'o'};
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		#this might glitch out on rare occasions...
		warn localtime(time())." [WARN] $0: no output vcf file specified or not -e 'out.vcf' so printing to STDOUT"."\n";
		$out =*STDOUT;
	}
	my $fieldstring = $opts -> {'f'};
	my $fields;
	@{$fields} = split(',', $fieldstring);
	for my $field (@{$fields}){
		my @fielddata = split('=',$field);
		undef($field);
		if(scalar(@fielddata) == 2){
			%{$field} = ('name' => $fielddata[0], 'oldname' => $fielddata[1]);
		}else{
			%{$field} = ('name' => $fielddata[0], 'oldname' => $fielddata[0]);
		}
		
	}
	#die Dumper($fields);
	warn localtime(time())." [INFO] $0: Parsing header.";
	my $vcf = Vcf->new(file=>$vcfin);
	#warn "parse header";
	$vcf->parse_header();
	my $vcfout = Vcf->new(file=>$vcfin);#to get header from input file
	$vcfout->parse_header();
	InfoHeaderToFormat('vcf' => $vcfout,'fields' => $fields,);
	$vcfout -> add_columns('FORMAT',$opts -> {'g'});
	#warn "get annotations";
	print {$out} $vcfout->format_header();

	#die "header check";
	#$vcf->recalc_ac_an(0);
	warn localtime(time())." [INFO] $0: Iterating records."; my $recordcount=0;
	while (my $x=$vcf->next_data_hash()){
#		if(defined($annheaderline) && scalar(@{$annheaderline})){
#			FillSNPEFFANNFields('record' => $x,'ann' => $annotations);
#		}
		InfoRecordToFormat('vcf' => $vcfout,'fields' => $fields,'genotype' => $opts -> {'g'}, 'record' => $x);
		#die Dumper $x;
		print {$out} $vcfout->format_line($x);
		#die Dumper $x;
		$recordcount++;
#               #if(scalar(@{GetAlt($x)})>1){
#		#		die Dumper($x)."\n".GetAlt($x)."\n".$. ;
#		#}

	}
	$vcf->close();
	warn localtime(time())." [INFO] $0: Done. Processed $recordcount";
}

sub InfoHeaderToFormat {
	my $self;
        %{$self}= @_;
        my $vcf = $self -> {'vcf'} or die "no vcf object as input";
	my $fields = $self -> {'fields'} or die "no caller string as input";
	for my $field (@{$fields}){
		#die Dumper($vcf);
		$vcf -> add_header_line({'key' => 'FORMAT',
			'ID' => $field -> {'name'},
			'Number' => $vcf -> {header} -> {INFO} -> {$field -> {'oldname'}} -> {'Number'},
			'Type' => $vcf -> {header} -> {INFO} -> {$field -> {'oldname'}} -> {'Type'} ,
			'Description' => $vcf -> {header} -> {INFO} -> {$field -> {'oldname'}} -> {'Description'}});
		$vcf -> remove_header_line(key=>'INFO', ID=>$field -> {'oldname'});
	}
}
sub InfoRecordToFormat {
	my $self;
        %{$self}= @_;
        my $vcf = $self -> {'vcf'} or die "no vcf object as input";
	my $fields = $self -> {'fields'} or die "no caller string as input";
	my $record= $self -> {'record'} or die "no record as input";
	my $genotype= $self -> {'genotype'} or die "no genotype(samplename) as input";
	for my $field (@{$fields}){
		if(defined($record -> {'INFO'} -> {$field -> {'oldname'}})){
			$vcf -> add_format_field($record,$field -> {'name'});
			$record -> {'gtypes'} -> {$self -> {'genotype'}} -> {$field -> {'name'}} = $record -> {'INFO'} -> {$field -> {'oldname'}};
			delete($record -> {'INFO'} -> {$field -> {'oldname'}});
			#die Dumper($record);
		}
	}
}
