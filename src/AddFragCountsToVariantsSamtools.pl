#!/usr/bin/env perl
use warnings;
use strict;
use Vcf;
use Getopt::Std;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
&main();
sub main {
	#use List::Util qw/max/; my $in;
	
	my $use="$0 -v vcf [-f 0.05 -s path/to/samtools] sample_1,fgbio.sample_1.FilteredConsensusVariants.bam [sample_n,fgbio.sample_n.FilteredConsensusVariants.bam] > annotated.vcf
	quick n dirtay fragmentcounts into your variants
";
	my $opts;
	$opts ={'f' => '0.08',#0.08 to be conservative 0.05 maybe better
		's' => 'samtools'};#idk if i need more para
	getopts('f:v:o:s:', \%{$opts});
	
	warn "## ". localtime(time()). " ## INFO ## Reading into mem the frag counts.\n";
	
	my $sampletobaminfo;
	for my $sampleandbaminfo (@ARGV){
		my ($samplename,$bam) = split(',',$sampleandbaminfo);
		$sampletobaminfo -> {$samplename} = $bam;
	}
	warn "## ". localtime(time()). " ## This was all that you can see".Dumper($sampletobaminfo)." ";
	
	#warn localtime(time())." [INFO] $0: Parsing header.";
	my $vcf = Vcf->new(file=>$opts -> {'v'}) 
		or die "## ". localtime(time()). " ## ERROR cannot open input vcf. ".$opts -> {'v'}." ";

	my $out;
	if(defined($opts -> {'o'})){
		my $vcfout = $opts -> {'o'};
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn localtime(time())." [WARN] $0: no output vcf file specified printing to STDOUT"."\n";
		$out =*STDOUT;
	}

	#warn "parse header";
	$vcf->parse_header();
	$vcf -> add_header_line({'key' => 'FORMAT', 'ID' => "FDP", 'Number' => 'A','Type' => 'String','Description' => "Fragment DP. Using the consensus read workflow these are the unique by start and end position consensus fragment depth supporting the alternate alleles."});
	$vcf -> add_header_line({'key' => 'FORMAT', 'ID' => "FDPF1R2", 'Number' => 'A','Type' => 'String','Description' => "FDP for read pair orientation F1R2."});
	$vcf -> add_header_line({'key' => 'FORMAT', 'ID' => "FDPF2R1", 'Number' => 'A','Type' => 'String','Description' => "FDP for read pair orientation F2R1."});
		
	print {$out} $vcf->format_header();
	warn localtime(time())." [INFO] $0: Iterating records."; my $records=0;
	my $metrics;
	while (my $x = $vcf->next_data_hash()){
		#if(defined($annheaderline) && scalar(@{$annheaderline})){
		#	FillSNPEFFANNFields('record' => $x,'ann' => $annotations);
		#}
		#warn "## ". localtime(time()). " ## Test record before. ".Dumper($x)." ";
		#if(defined($sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}})){	
		#warn $..join(';',keys(%{$x -> {'gtypes'}})) . Dumper($x);
		my $fixupflag = 0;
		for my $gt (keys(%{$x -> {'gtypes'}})){
			#warn $sampletobaminfo -> {$gt}.";".$x -> {'gtypes'} -> {$gt} -> {'AD'}." ";
			#samtools view  /scratch/umcg-mterpstra/projects/18FEB2021_editstrat1_healthycontrols_plasma/CallMolecularConsensusReads/104103-001-006.bam  -L <(echo -e '2\t61128285\t61128286') 			
			if(defined($sampletobaminfo -> {$gt}) && defined($x -> {'gtypes'} -> {$gt} -> {'AD'}) && not($x -> {'gtypes'} -> {$gt} -> {'AD'} =~ m!\.(/\.)+!) ){
				#warn $gt;
				my @advals = split(',',$x -> {'gtypes'} -> {$gt} -> {'AD'});
				#sum ad
				my $sum = 0;
				map{$sum +=$_;}(@advals);

				#choose to annotate if frequency > set value || 0.05 and if bigger then 1 for efficiency and reduction of garbage processing;
				my $annotate = 0;
				for my $ad (@advals[1..$#advals]){
					if($ad > 1 && $ad/$sum <= $opts -> {'f'}){
						$annotate++;
						$fixupflag++;
					}
					
				}
				if($annotate > 0){
					#here the gore starts					
					$metrics -> {'annotated genotypes'}++;
					$metrics -> {'pos'} -> { 
						$x -> {'CHROM'} .
						":".
						$x -> {'POS'} .
						"-".
						$x -> {'REF'} .
						">" .
						join(",", @{$x -> {'ALT'}})
						} ++;
					my $fragcounts = GetFragCounts($x,$sampletobaminfo -> {"$gt"});
					my (@fdp,@fdpF1R2,@fdpF2R1);
					for my $alt (@{$x -> {'ALT'}}){
						push(@fdp,$fragcounts -> {$alt} -> {'fdp'});
						push(@fdpF1R2,$fragcounts -> {$alt} -> {'fdpf1r2'});
						push(@fdpF2R1,$fragcounts -> {$alt} -> {'fdpf2r1'});
					}
					#die Dumper($x,$fragcounts) if(scalar(@fdp) > 10);
					$vcf->add_format_field($x,'FDP'); $x -> {gtypes} -> { $gt } -> {'FDP'}=join(',',@fdp);
					$vcf->add_format_field($x,'FDPF1R2'); $x -> {gtypes} -> { $gt } -> {'FDPF1R2'}=join(',',@fdpF1R2);
					$vcf->add_format_field($x,'FDPF2R1'); $x -> {gtypes} -> { $gt } -> {'FDPF2R1'}=join(',',@fdpF2R1);
				}
			}
		}
		if($fixupflag){
			for my $gt (keys(%{$x -> {'gtypes'}})){
				if(not(defined($x -> {gtypes} -> { $gt } -> {'FDP'}))){
					$x -> {gtypes} -> { $gt } -> {'FDP'} = '.';
					$x -> {gtypes} -> { $gt } -> {'FDPF1R2'} = '.';
					$x -> {gtypes} -> { $gt } -> {'FDPF2R1'} = '.';
				}
			}
		}
		
		#die "## ". localtime(time()). " ## Test record after. ".Dumper($x)." ";
		print {$out} $vcf->format_line($x);
		$records++;
		#;die if $. > 136;
	}
	$vcf->close();
	warn localtime(time())." [INFO] $0: Done. Processed $records";
}

################################################################################
#general stuff
sub CmdRunner {
	my $ret;
	my $cmd = join(" ",@_);
	
	warn localtime( time() ). " [INFO] system call:'". $cmd."'.\n";
	
	@{$ret} = `set -o pipefail && ($cmd )2>&1`;
	if ($? == -1) {
		die localtime( time() ). " [ERROR] failed to execute: $!\n";
	}elsif ($? & 127) {
		die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
		 ($? & 127),  ($? & 128) ? 'with' : 'without';
	}elsif ($? != 0) {
		die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
	         ($? & 127),  ($? & 128) ? 'with' : 'without';
	}else {
		warn localtime( time() ). " [INFO] " . sprintf "child exited with value %d\n", $? >> 8;
	}
	
	for (@{$ret}){
		$_ =~ s!^(\s*)!$1    $0:!g;
	}

	return $ret;
}
sub SamtoolsViewRunner {
	my $x = shift(@_);
	#my $start = shift(@_);
	#my $end = shift(@_);
	my $bam = shift(@_);
	#my $maxlength;map{$maxlength = lenght($_)if($maxlength < $_)}(@{$x -> {'ALT'}});
	my $bed = GetBedVariant($x);
	my $samlines = CmdRunner("echo -e '".$bed -> [0]."\\t".$bed -> [1]."\\t". $bed -> [2]."' | samtools view '$bam' -L /dev/stdin");
	#samtools view  /scratch/umcg-mterpstra/projects/18FEB2021_editstrat1_healthycontrols_plasma/CallMolecularConsensusReads/104103-001-006.bam  -L <(echo -e '2\t61128285\t61128286'
	return $samlines;
}
sub GetFragCounts {
	my $vcfrecord = shift @_;
	my $bam = shift @_;
	my $samreads = SamtoolsViewRunner($vcfrecord,$bam);
	#
	my $seqcounts;
	
	for my $samread (@{$samreads}){
		my $sam = GetSam($samread);
		my ($seq,$qual,@data) = GetReadBase($vcfrecord,$sam);
		$seqcounts -> {$seq} -> {'quals'} -> {$qual} ++; #<-magic here
		$seqcounts -> {$seq} -> {'counts'} ++; 
		$seqcounts -> {$seq} -> {'TLEN'} -> {GetSamFragmentInterval($sam)}++;
		#die Dumper($sam,GetSamFragmentInterval($sam)); 
	}
	#my @data =GetReadBase($vcfrecord,$samreads -> [0]);
	my $fragcounts;
	for my $alt (@{$vcfrecord -> {'ALT'}}){
		my ($fdp,$fdpf1r2,$fdpf2r1);		
		if(defined($seqcounts -> {$alt})){
			for my $tlen (keys(%{$seqcounts -> {$alt} -> {'TLEN'}})){
				$fdp -> {substr($tlen,0,length($tlen)-2)} ++;
				if(substr($tlen,length($tlen)-1,1) eq "+"){
					$fdpf1r2 -> {$tlen} ++;
				}else{
					$fdpf2r1 -> {$tlen} ++;
				}
			}
			
		}
		#warn Dumper($fdp,$fdpf1r2,$fdpf2r1)." ";
		$fragcounts -> {$alt} -> {'fdp'} = keys(%{$fdp});
		$fragcounts -> {$alt} -> {'fdpf1r2'} =keys(%{$fdpf1r2});
		$fragcounts -> {$alt} -> {'fdpf2r1'} =keys(%{$fdpf2r1});
	}
	return $fragcounts;
	#die Dumper($vcfrecord,$seqcounts,$fragcounts).$samreads -> [0],$samreads -> [-1],$samreads -> [-50] if (scalar());
	
}

################################################################################
#SAM Stuff needs to be checked
#
sub SamAsString {
	my $s = shift @_;
	my $string = join("\t",@{$s})."\n";
	return $string;
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
	#warn Dumper($_);
	for my $r (@{$_}){
		print join("\t",@{$r})."\n";
	}
}

sub CigarSamParser {
	my $s= shift(@_);
	my $cigar;
	my $c;
	@{$c} = split(/([A-Z=])/,GetSamCigarRead($s));
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
	my $c = shift(@_);
	#warn "CigarDump".Dumper(\$c)."";
	my $length = 0;
	while(my $operation = shift(@{$c})){
		my ($operator, $amount) = %{$operation};
		die "invalid thing".Dumper($operation) if(scalar(keys(%{$operation}))!=1);
		if($operator =~ /^[M\=X]$/){
			$length += $amount;
		}elsif($operator =~ /^[IS]$/){
			#nothing
		}elsif($operator =~ /^[DN]$/){
			$length += $amount;
		}elsif($operator =~ /^[HP]$/){
			die "Fatal Cigar operations \^[HP]\$ are not supported ".Dumper($c)." ";
		}
	}
	#warn "CigarDump".Dumper(\$c,\$length)."";
	return $length;
}
sub CigarParsedAsString {
	my $c = shift(@_);
	#warn "CigarDump".Dumper(\$c)."";
	my $string = '';
	while(my $operation = shift(@{$c})){
		my ($operator, $amount) = %{$operation};
		die "invalid thing".Dumper($operation) if(scalar(keys(%{$operation}))!=1);
		if($amount == 1 ){#best practice write as:
			$string .= $amount. $operator;
		}elsif($amount == 0 ){
			#ignore do not grow string
		}else{
			$string .= $amount. $operator;
		}
	}
	return $string;
}
sub GetSamNameRead {
	my $s= shift(@_);
	#die $s -> [0];
	return $s -> [0];
}

sub GetSamCigarRead {
	my $s= shift(@_);
	#confess " #############3here" if($s == 75);
	return $s -> [5];
}
sub GetSamFlagRead {
	my $s= shift(@_);
	#die $s -> [1];
	return $s -> [1];
}
sub GetSamChromRead {
	my $s= shift(@_);
	#confess " #############3here" if($s == 75);
	return $s -> [2];
}
sub GetSamPosRead {
	my $s= shift(@_);
	return $s -> [3];
}
sub SetSamPosRead {
	my $s= shift(@_);
	my $pos = shift(@_);
	$s -> [3] = $pos;
	return $s;
}
sub GetSamSeqRead {
	my $s= shift(@_);
	return $s -> [9];
}
sub SetSamSeqRead {
	my $s= shift(@_);
	my $seq = shift(@_);
	$s -> [9] = $seq;
	return $s;
}
sub GetSamQualRead {
	my $s= shift(@_);
	return $s -> [10];
}
sub GetSamMapqRead {
	my $s= shift(@_);
	return $s -> [4];
}

sub SetSamQualRead {
	my $s= shift(@_);
	my $qual = shift(@_);
	$s -> [10] = $qual;
	return $s;
}
sub GetSamRnextRead {
	my $s= shift(@_);
	return $s -> [6];
}
sub GetSamPnextRead {
	my $s= shift(@_);
	return $s -> [7];
}
sub GetSamTlenRead {
	my $s= shift(@_);
	return $s -> [8];
}

sub IsPrimaryAlignment {
	my $s=shift(@_);
	#warn "test".(GetFlagRead($r) & 256);
	if((GetSamFlagRead($s) & 256)){
		
		return 0;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 1;
	}
}
sub IsReverseAlignment {
	my $s=shift(@_);
	
	if((GetSamFlagRead($s) & 16)){
		
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 0;
	}
}

sub GetSamReadPairEnd {
	my $s=shift(@_);
	
	if((GetSamFlagRead($s) & 64 )){
		return 1;
	}elsif((GetSamFlagRead($s) & 128)){
		return 2;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return undef;
	}
}
sub GetSamOptionalFields{
	my $s=shift(@_);
	#my $s= shift(@_);
	#my $qual = shift(@_);
	#return $s -> [3];
	my $opt;
	#warn "$#{$s}";
	@{$opt} = @{$s}[11 .. -1];
	return $opt;
}
sub GetSamFragmentInterval {
	my $sam = shift(@_);
	my $interval;
	my ($chr,$start,$end,$orientation);
	$chr = GetSamChromRead($sam);
	#process paired primary alignment reads only for simplicity
	if(IsPrimaryAlignment($sam)){
		if(defined(GetSamReadPairEnd($sam))){
			if(IsReverseAlignment($sam)){
				if(GetSamReadPairEnd($sam) == 1){
					$start = GetSamPosRead($sam) + GetCigarLength(CigarSamParser($sam)) - abs(GetSamTlenRead($sam));
					$end = GetSamPosRead($sam) + GetCigarLength(CigarSamParser($sam));
					$orientation = "-";
				}else{
					$start = GetSamPosRead($sam) + GetCigarLength(CigarSamParser($sam)) - abs(GetSamTlenRead($sam));
					$end = GetSamPosRead($sam) + GetCigarLength(CigarSamParser($sam));
					$orientation = "+";
				}
			}else{
				if(GetSamReadPairEnd($sam) == 1){
					$start = GetSamPosRead($sam);#ok
					$end = GetSamPosRead($sam) + abs(GetSamTlenRead($sam));
					$orientation = "+";
				}else{
					$start = GetSamPosRead($sam);#ok
					$end = GetSamPosRead($sam) + abs(GetSamTlenRead($sam));
					$orientation = "-";
				}
			}
		}
	}
	$interval = $chr.":".$start."-".$end.":".$orientation;
	#.":Tlen".GetSamTlenRead($sam).":cig".GetSamCigarRead($sam).GetCigarLength(CigarSamParser($sam)).":pos".GetSamPosRead($sam).":readpair".GetSamReadPairEnd($sam).":Reverse".IsReverseAlignment($sam).":name".GetSamNameRead($sam);
	return $interval;
}
sub GetSam {
	#aka GetSam($samline));
	my $r = shift(@_);
	my $tsv;
	@{$tsv} = split("\t",$r);
	my $sam; 
	@{$sam}=(GetSamNameRead($tsv),		#0
		GetSamFlagRead($tsv),	#1
		GetSamChromRead($tsv),	#2
		GetSamPosRead($tsv),		#3
		GetSamMapqRead($tsv),	#4
		GetSamCigarRead($tsv),	#5
		GetSamRnextRead($tsv),	#6
		GetSamPnextRead($tsv),	#7
		GetSamTlenRead($tsv),	#8
		GetSamSeqRead($tsv),		#9
		GetSamQualRead($tsv),	#10
#		GetSamOptionalFields($r) 	#11+
	);
	push(@{$sam},GetSamOptionalFields($tsv)) if(defined(GetSamOptionalFields($tsv)));
	return $sam;
}
#sub CorrectIntervalForCigar {
#	
#}

sub GetReadBase {
	my $vcfrecord = shift @_;
	my $sam = shift @_;
	$vcfrecord -> {'POS'};
	my $sampos = GetSamPosRead($sam);
	my $cigar = CigarSamParser($sam);
	my $seq = GetSamSeqRead($sam),
	my $qual = GetSamQualRead($sam),
	
	my $loc = GetIntervalVariant($vcfrecord);
	$loc -> [1] = $loc -> [1] - $sampos;
	$loc -> [2] = $loc -> [2] - $sampos;
	#warn Dumper($loc)." ";
	my $testflanks = 0 ; #should be 0 for regular use 
	#my $testoffset = 0 ;
	my @startend;
	for my $pos (@{$loc}[1,2]){
		my $seqindex = 0;#where the start/end is resp
		my $cigarindex = 0;
		while($cigarindex < scalar(@{$cigar}) && $pos >= 0){
					
			my $ref= @{$cigar}[$cigarindex];
			#do some calculations and set SamPos,SamCigar...to the correct positionsand hard clip when needed.
			my ($operation,$amount) = %$ref;
			if($operation =~ /^[M\=X]$/){
				#matching alignment space
				if($pos-$amount>=0){
					#Incomplete trim
					$seqindex += $amount;
				}else{
					#complete trim
					$seqindex += ($pos);
				}
				
				$pos -= $amount;
			}elsif($operation =~ /^[IS]$/){
				#insertion in alignment space
				if($pos-$amount>=0){
					#Incomplete trim
					$seqindex += $amount;
				}else{
					#complete trim
					$seqindex += ($pos);
				}
			}elsif($operation =~ /^[DN]$/){
				#deletion in alignment space
				if($pos-$amount>=0){
					#Incomplete trim
				}else{
					#complete trim
				}
				$pos -= $amount;#so that we can remove the deletio
			}elsif($operation =~ /^[HP]$/){
				die "Fatal Cigar operations \^[HP]\$ are not supported ".Dumper($cigar,$sam,\$pos)." ";
			}
			die "Alignment position errror\n" if($loc < 0);
			#warn "### pos current:".$pos;
			$cigarindex++;
		}
		push(@startend,$seqindex);
	}
	#GetSamTlenRead($sam);
	return substr($seq,$startend[0] - $testflanks,$startend[1]-$startend[0] +1 + $testflanks +$testflanks),substr($qual,$startend[0] - $testflanks,$startend[1]-$startend[0] +1 + $testflanks +$testflanks),\@startend,$loc;
}
#sub GetSamMetrics {
#
#}

################################################################################
#vcf
sub GetIntervalVariant {
	my $vcfrecord = shift @_;
	my $interval;
	@{$interval} = ($vcfrecord -> {'CHROM'},($vcfrecord -> {'POS'}),($vcfrecord -> {'POS'} - 1 + length($vcfrecord -> {'REF'})));
	return  $interval;
};
sub GetBedVariant {
	my $vcfrecord = shift @_;
	my $bed;
	$bed = GetIntervalVariant($vcfrecord);
	$bed -> [1] = $bed -> [1] - 1;#format conversion
	return  $bed;
};

