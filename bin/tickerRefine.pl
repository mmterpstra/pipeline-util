#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

main();

sub main{
	#multihits on probes something needs to care of this. This is it.
	#use:
	#paste \
	# <(bedtools intersect -wao -bed -a <(\
	#  samtools view -Sb 3294553_IN-14_1545_A1_1_simple.sam) \
	#  -b /groups/umcg-oncogenetics/tmp04/res/probe_ET1262F_2_182.ensembl.bed -S |\
	# perl tickerRefine.pl - ) <(grep -vP '^\@'  3294553_IN-14_1545_A1_1_simple.sam  ) |\
	#tail -n 40
	
	#tail -n 40 || tickertape.pl
	
	warn "Commandline: $0 ".join(' ',@ARGV)."\n";
	
	my $recordLast;
	my $buffer;
	my $itstimetodie=5;
	
	while(<>){
		my $record;
		$record = DumbReader($_);

		if(defined($recordLast) && EqualReadRecords($recordLast,$record)){

			#buffer some
			push(@{$buffer}, $recordLast);

		}elsif(defined($buffer) && not(EqualReadRecords($recordLast,$record))){

			push(@{$buffer}, $recordLast);
			
			#warn Dumper($buffer)."#";
			
			$recordLast = RefineOverlapsBuffer($buffer);
			
			#die Dumper($buffer,$record,$recordLast)if($itstimetodie<0);
			#$itstimetodie--;
			
			print Writer($recordLast);
			
			$buffer=undef;
						
		
		#also read not(EqualReadRecords($recordLast,$record)) below, but not inserted because efficiency
		}elsif(not(defined($buffer)) && defined($recordLast)){

			#push(@{$buffer}, $recordLast);
			print Writer($recordLast);
			
		}elsif(not(defined($buffer)) && not(defined($recordLast))){
			#just exit and assign $recordLast = $record;
		}else{

			die 'Invalid something';

		}
		
		$recordLast = $record;
	}
	#purge buffer
	print Writer($recordLast);
}

sub EqualReadRecords{
	my $r=shift(@_);
	my $r2=shift(@_);
	if((GetChrRead($r) eq GetChrRead($r2)) 
		&& (GetStartRead($r) == GetStartRead($r2)) 
		&& (GetEndRead($r) == GetEndRead($r2)) 
		&& (GetNameRead($r) eq GetNameRead($r2)) 
		&& (GetStrandRead($r) eq GetStrandRead($r2))){
		return 1;
	}
	return 0;
}
sub GetChrRead{
	my $r=shift(@_);
	return($r -> [0]);
}
sub GetStartRead{
	my $r=shift(@_);
	return($r -> [1]);
}
sub GetEndRead{
	my $r=shift(@_);
	return($r -> [2]);
}
sub GetNameRead{
	my $r=shift(@_);
	return($r -> [3]);	
}
sub GetStrandRead{
	my $r=shift(@_);
	return($r -> [5]);	
}

sub GetChrProbe{
	my $r=shift(@_);
	return($r -> [12]);
}
sub GetStartProbe{
	my $r=shift(@_);
	return($r -> [13]);
}
sub GetEndProbe{
	my $r=shift(@_);
	return($r -> [14]);
}
sub GetNameProbe{
	my $r=shift(@_);
	return($r -> [15]);	
}
sub GetStrandProbe{
	my $r=shift(@_);
	return($r -> [17]);	
}

sub DumbReader{
	$_= shift(@_);
	chomp;
	my $record;
	@{$record}=split("\t");
	return $record;
}
sub Writer{
	$_= shift(@_);
	join("\t",@{$_})."\n";
}

sub GetHeader{
	my $record = shift(@_);
	#warn $record -> [3];
	return($record -> [3]);# or die 'Record does not contain this many fields!'.Dumper($record);
}
sub GetH2{
	my $record = shift(@_);
	#warn $record -> [19];
	if($record -> [0] ne "."){
		return($record -> [19]);# or die 'Record does not contain this many fields!';
	}else{
		#warn $record -> [20];
		return($record -> [20]);# or die 'Record does not contain this many fields!';
	}
}

sub RefineOverlapsBuffer {
	my $b = shift(@_);
	
	my $r1= shift(@{$b});
	if(scalar(@{$b}) > 0){
		while(my $r2= shift(@{$b})){
			$r1=GetBestOverlap($r1,$r2);
		}
	}
	push(@{$b},$r1);
	return $r1;
}
sub GetBestOverlap{
	my $r1= shift(@_);
	my $r2= shift(@_);
	
	#has two record or more:
	if(Get3PrimeOverlap($r1)<Get3PrimeOverlap($r2)){
		$r1=$r2;
	}
	#die Dumper($r1,$r2)."#";
	return $r1;
}

#notes for sub below
#'12',
#          '56493411',
#          '56493512',
#          'HISEQ-MFG:688:C7Y7WACXX_TCCAAT:1:1101:12442:2813',
#          '60',
#          '-',


#'12',
#          '56493381',
#          '56493431',
#          'ERBB3|NM_001982_exon_23_0_chr12_56493432_f_0_1937047',
#          '1',
#          '+',

#'12',
#          '56493507',
#          '56493557',
#          'ERBB3|NM_001982_exon_24_0_chr12_56493622_f_0_2022598',


#56493381	56493411			56493431	56493507	56493512	56493557
#|			|					|			|			|			|
#>>>>>>>>>>>>>>p1>>>>>>>>>>>>>>>>
#			<<<<<<<<<<<<<<<<<<<<<<read<<<<<<<<<<<<<<<<<<<
#											>>>>>>>>>>>>>>p2>>>>>>>>>
#
##has any 3p overlap -
#pend >= readstart + 1 
#pstart - wiggle <= readstart
#
#aka
#GetEndProbe($r) >= GetStartRead($r) + 1;
#GetStartProbe($r) - $wiggle <= GetStartRead($r);
#$overlap=GetEndProbe($r) - GetStartRead($r) + 1;
#56493381	56493411			56493431	56493507	56493512	56493557
#|			|					|			|			|			|
#<<<<<<<<<<<<<<p1<<<<<<<<<<<<<<<<
#			>>>>>>>>>>>>>>>>>>>>>>read>>>>>>>>>>>>>>>>>>>
#											<<<<<<<<<<<<<<p2<<<<<<<<<
#
##has any 3p overlap +
#pend + wiggle >= readend
#pstart + 1<= readend
#
#aka
#GetEndProbe($r) + $wiggle >= GetEndRead($r);
#GetStartProbe($r) +1 <= GetEndRead($r);
#$overlap=GetEndRead($r)-GetStartProbe($r) + 1;

sub Get3PrimeOverlap{
	my $r= shift(@_);
	my $overlap=0;
	my $wiggle = 0;
	
	if(	GetStrandRead($r) ne GetStrandProbe($r) 
		&& GetStrandRead($r) ne '.'
		&& GetStrandProbe($r) ne '.'){
		
		if(GetStrandRead($r) eq '+'
			&& GetEndProbe($r) + $wiggle >= GetEndRead($r)
			&& GetStartProbe($r) +1 <= GetEndRead($r)){
				
			#--->read>
			#-----<probe<
			
			$overlap = GetEndRead($r) - GetStartProbe($r);
			#die Dumper(\$overlap,\$overlap,$r)."#";
		}elsif(GetStrandRead($r) eq '-'
			&& GetEndProbe($r) >= GetStartRead($r) + 1 
			&& GetStartProbe($r) - $wiggle <= GetStartRead($r)){
				
			#---<read<
			#->probe>
			
			$overlap = GetEndProbe($r) - GetStartRead($r);
			#die Dumper(\$overlap,\$overlap,$r);
		}
	}else{
		#warn Dumper($r, \$overlap);
		#die Dumper($r, \$overlap)."#" if($overlap < 0);
	}
	#warn Dumper($r, \$overlap)."#";
	#die Dumper($r, \$overlap)if($overlap > 0);
	
	return $overlap;
}


