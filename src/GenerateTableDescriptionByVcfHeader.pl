#!/usr/bin/env perl -w

use warnings;
use strict;
use Vcf;
use Data::Dumper;
my $use =<<"END";

$0 file.vcf > out.table
Generates tabular description of vcf content

END
#

die "Invalid vcf file specified!! '$ARGV[0]'\n $use\n$!" if(not(-e $ARGV[0]));
VcfHeaderToTable($ARGV[0]);

sub VcfHeaderToTable{
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	#if($_[1] && not(-e $_[1])){
	#	my $vcfout = $_[1];
	#	open($out,'>',$vcfout) or die 'Cannot write vcf file';
	#}else{
	#	warn 'warning: no out.vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	#}
	
	
	my $vcf = Vcf->new(file=>$vcfin) or die "cannot open vcf file $vcfin\n";
	$vcf->parse_header();
	#print Dumper($vcf)."\n";
	print "Description and Filtering results of file:'$vcfin'\n";
	print "The following fields are present having info about the variant:"."\n";
		print "\t".'ID'."\t".'Type'."\t".'Description'."\n";
	
	for my $headerline (@{$vcf->{'header_lines'}}){
		print "\t".$headerline->{'ID'}."\t".$headerline->{'Type'}."\t".$headerline->{'Description'}."\n"if($headerline->{'key'} eq 'INFO') ;
	}
	print "The following fields are present having info about the genotype:"."\n";
			print "\t".'ID'."\t".'Type'."\t".'Description'."\n";
	for my $headerline (@{$vcf->{'header_lines'}}){
		print "\t".$headerline->{'ID'}."\t".$headerline->{'Type'}."\t".$headerline->{'Description'}."\n"if($headerline->{'key'} eq 'FORMAT') ;
	}
	print "The following Filters were applied:"."\n";
	print "\t".'ID'."\t".'Description'."\n";
	for my $headerline (@{$vcf->{'header_lines'}}){
		print "\t".$headerline->{'ID'}."\t".$headerline->{'Description'}."\n"if($headerline->{'key'} eq 'FILTER') ;
	}
	
	#print {$out} $vcf->format_header();
	$vcf->recalc_ac_an(0);
	my %filterdata;
	
	$filterdata{'Total'}=0;

	while (my $x=$vcf->next_data_hash()){
		$filterdata{'Total'}++;
		for my $filter (@{$x->{'FILTER'}}){
			$filterdata{'SingleCounts'}{$filter}++;
		}
		$filterdata{'Complextable'}{join(',',(@{$x->{'FILTER'}}))}++;
		#for my $F (keys(%{$x->{'INFO'}})){
		#	next if($F =~ /_[0123456789]$/);
		#}
		#for my $GF (@{$x->{'FORMAT'}}){
		#}
		#print {$out} $vcf->format_line($x);
	#FILTER	
		
#		 if ($. > 2);
		
	}
	$vcf->close();
	
	print "Count by single filter:"."\n";
	print "\t".'Filtername'."\t".'Count'."\n";
	
	print "\t".'Total'."\t".$filterdata{'Total'}."\n";
	if(defined($filterdata{'SingleCounts'}) && defined(getHashRef($filterdata{'SingleCounts'}))){
		for my $filterName (keys(%{getHashRef($filterdata{'SingleCounts'})})){
			print "\t".$filterName."\t".$filterdata{'SingleCounts'}{$filterName}."\n";
		}

		print "Count by combined filters:"."\n";
		print "\t".'Filtername(Comma,Separated)'."\t".'Count'."\n";
		print "\t".'Total'."\t".$filterdata{'Total'}."\n";
	
		for my $filterName (keys(%{getHashRef($filterdata{'Complextable'})})){
			print "\t".$filterName."\t".$filterdata{'Complextable'}{$filterName}."\n";
		}
	}else{
		print "No filter data present pls populate the filter field or the variants.";
	}
}
sub getHashRef{
	my $href = $_[0];
	return $href;
}

