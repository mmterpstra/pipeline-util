#!/usr/bin/env perl
use warnings;
use strict;
use Vcf;
#use Devel::Confess;
use Data::Dumper;
use pipeline::util;
use Getopt::Std;
my %opts;
getopt('cfn', \%opts);

if(not(defined($opts{'c'}))){
	$opts{'c'}=4;
}
if(not(defined($opts{'f'}))){
	$opts{'f'}=0.05;
}
if(not(defined($opts{'n'}))){
	$opts{'n'}=0.4;
}
if(not(defined($opts{'m'}))){
	$opts{'m'}=2;
}
my $use = <<"END";
$0 [options] TARGETVCF [VCF]
	Filters Target vcf with one or more VCF files from input
	-c int		Takes integer arg for min count difference between other samples
	-f float	Takes a float for min allele difference between other samples
	-n float	Max frequency in normal [defaults to 0.4]
	-m float	Frequency multiplier compared to normal (ftumor > fnormal * multiplier) [defaults to 2]
END

my $targetvcf = shift @ARGV or die localtime(time())." [ERROR] $0: No TARGETVCF option supplied. ARGV=[@ARGV]\n$use";

die "no valid 'VCF' specified on command line\n$use VCF=[@ARGV]" if( scalar(@ARGV)==0 || not( -e $ARGV[0]));

main("target"=>$targetvcf, "vcfs"=>\@ARGV ,'opts'=> \%opts);

sub main{
	my $self;
	%{$self} = @_;
	my $vcfin =$self -> {'target'} or die "$use";
	die "TARGETVCF does not exist. Arg:'$vcfin' " if(not(-e $vcfin));
	my $out =*STDOUT;

	warn localtime(time())." [INFO] $0: Parsing header of 'vcfin'";
	
	my $resource;
	die localtime(time())." [ERROR] $0: No resource vcf given.\n$use" if (scalar(@ARGV)== 0 );
	@{$resource -> {'invcfs'}} = @{$self -> {"vcfs"}};
	warn Dumper({'targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}}});
	my $walkdata;
	$walkdata = NewWalk('targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}});
	
	warn localtime(time())." [INFO] $0: Creating header of output";
		
	my $headerlines;
	@{$headerlines} = ({key=>'FILTER', ID=>'NoAD', Description=>'No Allele Depth values to filter on'}, 
		{key=>'FILTER', ID=>'ADDeltaCountlt'.$opts{'c'}, Description=>'LowADCount'},
		{key=>'FILTER', ID=>'ADDeltaFrequencylt'.$opts{'f'}, Description=>'LowADfrequency'},
		{key=>'FILTER', ID=>'ADNormFrequency'.$opts{'n'}, Description=>'Too high frequency in normal indicating probably het or hom in normal samples'},
		{key=>'FILTER', ID=>'ADMultFrequency'.$opts{'m'}, Description=>'Frequency multiplier compared to normal (ftumor < fnormal * multiplier) [defaults to '.$opts{'m'}.']'},
		{key=>'FILTER', ID=>'PASS', Description=>'All filters passed'});
	print FormatWalkTargetLineAsVcfHeader('walk' =>  $walkdata,'headerlines' => $headerlines);
	while(WalkToNext('walk'=> $walkdata)){
		#die Dumper($walkdata)."here"
		$walkdata = ADFilterTargetRecords('walk'=> $walkdata,'deltafrequency'=> $opts{'f'},'deltacount'=> $opts{'c'},'normfreq'=> $opts{'n'},'normfreqmult'=> $opts{'m'});
		print FormatWalkTargetLineAsVcfLine('walk'=> $walkdata);
	}
	#print "Last\n";
	$walkdata = ADFilterTargetRecords('walk'=> $walkdata,'deltafrequency'=> $opts{'f'},'deltacount'=> $opts{'c'},'normfreq'=> $opts{'n'},'normfreqmult'=> $opts{'m'});
	print FormatWalkTargetLineAsVcfLine('walk'=> $walkdata);
	warn _formatwalkasvcflineswithfile('walk' => $walkdata);
	#TargetVcfReAnnotator('targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}});
}
