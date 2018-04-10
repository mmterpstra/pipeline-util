#!/usr/bin/env perl -w
use warnings;
use strict;
use Vcf;
use List::Util qw(max min);
use List::MoreUtils qw(uniq);
use Data::Dumper;

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
	}elsif(-e $_[1]){

		my $vcfout = $_[1];
		warn "## Overwriting outfile '".$vcfout."'";
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn 'warning: no output vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	}
	
	
	my $vcf = Vcf->new(file=>$vcfin);
	$vcf->parse_header();
	print {$out} $vcf->format_header();
	#$vcf->recalc_ac_an(0);
	while (my $x=$vcf->next_data_hash()){
		#warn $x->{'INFO'}{'VariantType'}; 
		my @samples = GetSamples($x);
		#warn Dumper(@samples); 
		for my $sample (@samples){
			#get funky record based on the fact that i seen them as GT = 1/1 and a not-compliant AD
			if((HasFunkyGTPL('record' => $x, 'sample' => $sample) == 1)){
				#
				warn Dumper($x)." ";
				DeleteSampleFormatFields('record' => $x, 'sample' => $sample);
				#die Dumper($x)." ";
			}
			
		}
		
		print {$out} $vcf->format_line($x);
		#die;
		#die if ($. > 50);
	}
	$vcf->close();
}
sub GetSamples {
	my $self = shift @_ or die "No input record";
	return keys(%{$self->{gtypes}});
}
sub GetAD {
	my %self = @_;
	#
	return $self{'record'}->{'gtypes'} -> { $self{'sample'} } -> {'AD'};
	
}
sub GetPL {
	my %self = @_;
	#
	return $self{'record'}->{'gtypes'} -> { $self{'sample'} } -> {'PL'};
	
}
sub GetGT {
	my %self = @_;
	return $self{'record'}->{'gtypes'} -> { $self{'sample'} } -> {'GT'};
}
sub HasFunkyGTPL {
	my %self = @_;
	my $x = $self{'record'};
	my $sample = $self{'sample'};
	warn $sample;
	# by reverse reasoning it it mutch faster to find broken GT when it only happens in the 1/1 tag. 
	if(GetAD('record' => $x, 'sample' => $sample) &&
		not(GetAD('record' => $x, 'sample' => $sample) =~ m/^\.,\.$/)){
		my @AD = split(',',GetAD('record' => $x, 'sample' => $sample));
		my @ADsort = sort { $b <=> $a }(@AD);
		#warn Dumper(\@AD,\@ADsort);
		if(scalar(@AD) > 3 && GetSet('record' => $x) eq 'freebayes' && (GetGT('record' => $x, 'sample' => $sample) eq '1/1') ){
			if( (max(@AD) == $AD[1]) && ($ADsort[0] != $ADsort[1])){
				return 0;
			}else{
				return 1;
			}
		}
	}
	#if(GetPL('record' => $x, 'sample' => $sample) ){
	#	my @PL = split(',',GetPL('record' => $x, 'sample' => $sample));
	#	if(GetSet('record' => $x) eq 'freebayes' && (GetGT('record' => $x, 'sample' => $sample) eq '1/1') ){
	#		if( (min(@PL)==0) && ($PL[2] == 0)){
	#			return 0;
	#		}else{
	#			return 1;
	#		}
	#	}
	#}
	return 0;
}
sub DeleteSampleFormatFields {
	my %self = @_;
	$self{'record'}->{'gtypes'} -> { $self{'sample'} } = {'GT' => './.'};
	return $self{'record'};
}
sub GetSet {
	my %self = @_;
	return $self{'record'}->{'INFO'} -> {'set'};
}
#('record' => $x, 'sample' => $sample);
