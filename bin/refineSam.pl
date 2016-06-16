#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use Carp;

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
	
	my $recordsLast;
	my $buffer;
	
	while(<>){
		next if(m/^\@/);
		
		my $record;
		$record = DumbReader($_);

		if(defined($recordsLast) && EqualReadNames($recordsLast,$record)){

			#buffer some
			push(@{$buffer}, @{$recordsLast});

		}elsif(defined($buffer) && not(EqualReadNames($recordsLast,$record))){

			push(@{$buffer}, @{$recordsLast});
			
			#warn Dumper($buffer)."#";
			
			#die Dumper($buffer).$. ;#if(scalar(@{$buffer}) > 2);
			$recordsLast = PickBest($buffer);
			
			#die Dumper($buffer,$record,$recordLast)if($itstimetodie<0);
			#$itstimetodie--;
			
			print Writer($recordsLast);
			$recordsLast=undef;
			
			$buffer=undef;
						
		
		#also read not(EqualReadRecords($recordLast,$record)) below, but not inserted because efficiency
		}elsif(not(defined($buffer)) && defined($recordsLast)){

			#push(@{$buffer}, $recordLast);
			print Writer($recordsLast);
			$recordsLast=undef;
		}elsif(not(defined($buffer)) && not(defined($recordsLast))){
			#just exit and assign $recordLast = $record;
		}else{

			die 'Invalid something';

		}
		
		$recordsLast -> [0] = $record;
	}
	#purge buffer
	push(@{$buffer}, @{$recordsLast});
	$recordsLast = PickBest($buffer);
	print Writer($recordsLast);
}


sub main2{
	#multihits on probes something needs to care of this. This is it.
	#use:
	#paste \
	# <(bedtools intersect -wao -bed -a <(\
	#  samtools view -Sb 3294553_IN-14_1545_A1_1_simple.sam) \
	#  -b /groups/umcg-oncogenetics/tmp04/res/probe_ET1262F_2_182.ensembl.bed -S |\
	# perl tickerRefine.pl - ) <(perl refineSam.pl  3294553_IN-14_1545_A1_1_simple.sam  ) |\
	#tail -n 40
	#tail -n 40 || tickertape.pl
	
	warn "Commandline: $0 ".join(' ',@ARGV)."\n";
	
	my $recordLast;
	my $buffer;
	
	while(<>){
		#process/skip header	
		#print if(m/^\@/);
                next if(m/^\@/);

		#process reads
		my $record;
		$record = DumbReader($_);
		
		#warn Dumper($record).Writer($record);
		#die Dumper((GetChrRead($record),
		#	 GetStartRead($record),
		#	 GetEndRead($record), 
		#	 GetNameRead($record), 
		#	 GetStrandRead($record)));
		#die if($.== 200000);
	
		if(defined($recordLast) && EqualReadsAndMappings($recordLast,$record)){

			##buffer some
            push(@{$buffer}, $recordLast);
			#warn 'lastr:'.$recordLast.';r:'.$record;
			#warn Dumper($recordLast),Dumper($record).'here ';
			
		}elsif(defined($buffer) && not(EqualReadsAndMappings($recordLast,$record))){

			push(@{$buffer}, $recordLast);
			
			#warn Dumper($buffer)."#";
			
			$recordLast = PickBest($buffer);
			
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

			die 'Invalid something';#never happens

		}
		
		$recordLast = $record;
	}
	#purge buffered record
	print Writer($recordLast);
}
sub PickBest{
	my $b = shift(@_);
	
	if(HasPairedEndTag($b -> [0])){
		#run for /2 then /1
		my $r1;
		my $r2;
		if($r2=ShiftBufferReadEnd2($b)){
			if(scalar(@{$b}) > 0){
				while(my $r2next = ShiftBufferReadEnd2($b)){
					$r2=GetBestOverlap($r2,$r2next);
				}
			}	
		}
		if($r1=ShiftBufferReadEnd1($b)){
			my $r1old = $r1;
			if(scalar(@{$b}) > 0){
				while(my $r1next = ShiftBufferReadEnd1($b)){
					$r1=GetBestOverlap($r1,$r1next);
				}
			}
		}
		my $bNew;
		@{$bNew} = ($r1,$r2);
		
		return($bNew);
		
	}else{
		my $r1= shift(@{$b});
		
		if(scalar(@{$b}) > 0){
			while(my $r1next= shift(@{$b})){
				$r1=GetBestOverlap($r1,$r1next);
			}
		}
		push(@{$b},$r1);
		return @{$r1};
	}
	
}

sub EqualReadNames {
	my $r=shift(@_);
	if(scalar(@{$r->[0]})){
		#warn "scalar ref reassign of because of nesting ".Dumper($r,$r2);
		$r = $r -> [0];
		#warn "after:". Dumper($r,$r2);
		#die;
	};
	my $r2=shift(@_);
	if((GetNameRead($r) eq GetNameRead($r2)) ){
		return 1;
	}
	return 0;
	
}
sub EqualReadsAndMappings{
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
	return($r -> [2]);
}
sub GetStartRead{
	my $r=shift(@_);
	return($r -> [3]);
}
sub GetEndRead{
	my $r=shift(@_);
	my $rEnd;
	if(GetCigarRead($r) ne '*'){
		$rEnd = GetStartRead($r) + GetCigarLength($r) - 1;
	}else{
		$rEnd = 0;
	}
	#die $rEnd;
	return($rEnd);
}
sub GetNameRead{
	my $r=shift(@_);
	return($r -> [0]);	
}

sub GetStrandRead{
	my $r=shift(@_);
	my $ret;
	if(IsReverseAlignment($r)){
		$ret='-';
	}else{
		$ret='+';
	}
	return($ret);	
}
sub IsReverseAlignment {
	my $r=shift(@_);
	
	if((GetFlagRead($r) & 16)){
		
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 0;
	}
}

sub HasPairedEndTag {
	my $r=shift(@_);
	
	if((GetFlagRead($r) & 1)){
		
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 0;
	}

}
sub GetFlagRead{
	my $r=shift(@_);
	my $ret;
	$ret = $r -> [1];# or die 'Record does not contain this many fields!';
	#
	defined($ret) or confess "Invalid record at line ". $. .": ".Dumper($r);
	return($ret);
	#return($r -> [20]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetCigarRead{
	my $r=shift(@_);
	my $ret;
	$ret = $r -> [5] or die "Invalid record at line ". $. .": ".Dumper($r);
	#
	return($ret);
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
	confess Dumper($_) if(not(defined($_ -> [0])));
	for my $r (@{$_}){
		print join("\t",@{$r})."\n";
	}
}
sub IsPrimaryAlignment {
	my $r=shift(@_);
	#warn "test".(GetFlagRead($r) & 256);
	if((GetFlagRead($r) & 256)){
		
		return 0;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 1;
	}
}


sub GetBestOverlap{
	my $r1= shift(@_);
	my $r2= shift(@_);
	
	#has two record or more:
	if(GetFlagRead($r1)>GetFlagRead($r2)){
		$r1=$r2;
	}
	#die Dumper($r1,$r2)."#";
	return $r1;
}

sub CigarReader {
	my $r= shift(@_);
	my $c;
	my $cigar;
	@{$c} = split(/([A-Z=])/,GetCigarRead($r));
	while(my $tag = shift(@{$c})){
		if(looks_like_number($tag)){
			my %h = (shift(@{$c}) => $tag);
			push(@{$cigar}, \%h);
		}else{
			my %h = ($tag => 1);
			push(@{$cigar}, \%h);
		}
	}
	return $cigar;
}
sub GetCigarLength {
	my $r= shift(@_);
	my $clen = 0;
	my $cigar= CigarReader($r);
	for my $tag (@{$cigar}){
		#die Dumper($tag);
		my ($operation)= keys(%{$tag});
		my ($amount)= values(%{$tag});
		#warn "op:'$operation';amnt:'$amount';";
		if($operation =~ m/^[MN\=X]$/){
			$clen += $amount;
			
		}elsif($operation =~ m/^[IS]$/){
			#nothing
		}elsif($operation =~ m/^[D]$/){
			#only reduce overlap not increment trimbases left
			$clen += $amount;
		}else{
			die "invalid cigar operator '$operation'" . Dumper($r).Dumper($tag);
		}
	}
	#warn $clen ."#".Dumper($cigar)."here";
	
	return $clen;	
}

sub ShiftBufferReadEnd {
	my $b = shift(@_);
	#
	my $readEnd = shift(@_);#1 or 2
	#die "Funky data yes\n".Dumper($b) if($readEnd != 1 || $readEnd != 2 );
	#
	#warn "buffer contains:".Dumper $b;
	for(my $bIndex = 0; $bIndex < scalar(@{$b}); $bIndex++){
	#	#die "Funky data yes\n".Dumper($b) if(GetReadEnd($b -> [$bIndex]) != 1 || GetReadEnd($b -> [$bIndex]) != 2 );
	#	#return(splice(@{$b},$bIndex,1))if(GetReadEnd($b -> [$bIndex]) eq $readEnd);
		#my @tmpBuffer = shift(@{$b});
		if(GetReadEnd($b -> [$bIndex]) eq $readEnd){
	#		warn "buffer return contains: \n".Dumper($b).Dumper($b -> [$bIndex],$.);
			#too experimental... goes with return			
			return(splice(@{$b},$bIndex,1));
		}
	}
	
	return(undef);
	
}
sub ShiftBufferReadEnd1 {
	return(ShiftBufferReadEnd($_[0],1));
}
sub ShiftBufferReadEnd2 {
	return(ShiftBufferReadEnd($_[0],2));
}

sub GetReadEnd {
	my $r= shift(@_);
	my $readend;
	
	if((GetFlagRead($r) & 64)){
		
		return 1;
	}elsif((GetFlagRead($r) & 128)){
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 2;
	}else{
		confess "Invalid Pe Tag in flag! Flag=".GetFlagRead($r). Dumper($r)."line".$.;
	}
}

