#!/usr/bin/env perl -w
use warnings;
use strict;
use Vcf;
use Data::Dumper;
# input test
#$ARGV[0]='/home/terpstramm/Downloads/s12-193.annotated.vcf';
# use
my $use = <<"END";
	$0 in.vcf out.vcf
	ads Z scores based on the AD fields per sample
END
die "no valid 'in.vcf' specified on command line\n$use" if(not(defined($ARGV[0])) ||not -e $ARGV[0]);
#main
main($ARGV[0],$ARGV[1]);

sub main{
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
	print {$out} $vcf->format_header();
	#$vcf->recalc_ac_an(0);
	my $qssstats;
	while (my $x=$vcf->next_data_hash()){
		#warn $x->{'INFO'}{'VariantType'};
		#if(not($x->{'INFO'}{'VariantType'} =~ m/^MULTIALLELIC/)){
		my @samples = keys(%{$x->{gtypes}});
		warn Dumper($x) if ($x -> {'POS'} eq '151789588' ); 
		for my $sample (@samples){
				if(defined($x->{'gtypes'} -> {$sample} -> {'QSS'})){
					my @qss = split(',', $x -> {'gtypes'} -> {$sample} -> {'QSS'});
					$qssstats -> {scalar@{$x -> {ALT}}} -> {scalar(@qss)} ++;
					if(scalar(@{$x -> {ALT}}) +1 > scalar(@qss)){
						warn 'mismatching qss field found at '.$x -> {'CHROM'}.':'.$x -> {'POS'};
						while((scalar(@{$x -> {ALT}}) ) > scalar(@qss)){
							push(@qss,0);
						}
						$x->{'gtypes'} -> {$sample} -> {'QSS'} = join(',',@qss);
						$x->{'gtypes'} -> {$sample} -> {'MuTect2QSS'} =  join(',',@qss);
					}
				}
		}
		warn Dumper($x) if ($x -> {'POS'} eq '151789588' );
		#}else{
		#	warn Dumper($x);die;
		#}
		
		#warn Dumper($vcf);
		#warn Dumper($x);
		
		print {$out} $vcf->format_line($x);
		#die;
		#die if ($. > 50);
	}
	warn Dumper($qssstats);
	$vcf->close();
}
