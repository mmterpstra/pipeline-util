#!/usr/bin/env perl
use warnings;
use strict;
use Vcf;
#use Devel::Confess;
use Data::Dumper;
use pipeline::util;

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
	warn Dumper({'targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}}});
	my $walkdata;
	$walkdata = NewWalk('targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}});
	while(WalkToNext('walk'=> $walkdata)){
		Dumper($walkdata);
	}
	print "Last\n";
	_formatwalkasvcflineswithfile('walk' => $walkdata);
	#TargetVcfReAnnotator('targetvcf'=> $vcfin,'vcflist'=> \@{$resource -> {'invcfs'}});
}
