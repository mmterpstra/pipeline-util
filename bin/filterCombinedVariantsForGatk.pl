#!/usr/bin/perl
use warnings;
use strict;
my $last = '1'; 
my $last4;
my $_last;
while(<>){
	my @F = split("\t",$_);
	if(defined($F[1]) && $F[1] eq $last ){
		if($_ =~ m/set=GATK/){
			$_last =$_;
			$_ = '';
		}else{
			$_ = '';
		}
	};
	print $_last;
	$last=$F[1]; $_last=$_;$last4=$F[4];
}
print $_;
