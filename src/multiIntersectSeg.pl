#!/usr/bin/env perl -w
use warnings;
use strict;
use Data::Dumper;
#convert to bed
main();

sub main {
	my $srcdir = GetSourceDir();
	
	
	warn "use: perl $0 genome.dict /path/to/outdir files.seg \nMerges different seg files into single table. Juste like `bedtools multiinter` wich is used under the hood.";
	
	my @bedFiles=();
	
	my $dict = shift @ARGV;
	my $outDir = shift @ARGV;
	
	my $chrLenDat=ReadDict($dict);
	#die Dumper($chrLenDat);
	#my $offset=$chrLenDat->{'ChrOff'}->{$chrom}
	for my $segFile (@ARGV){
		my $seg=SegNew('file'=>$segFile);
		my $bedFile=$segFile;
		
		if(not($bedFile=~s/\.seg$/\.bed/)){
			$bedFile.=".bed";
		}
		
		warn localtime(time())." [INFO] Converting $segFile into bed\n";
		
		open(my $bedHandle,'>',$bedFile) or die "cannot write bed file".$bedFile;
		
		while($seg=SegNextDataHash($seg)){
			#warn Dumper($seg);
			print $bedHandle SegDataHashAsBedString($seg);
		}
		close($bedHandle);
		push(@bedFiles,$bedFile);
	}
	warn Dumper(\@bedFiles);
	my $dataMultiInter = MultiInter(\@bedFiles);
	
	#print string/table $dataMultiInter;
	open(my $mergedHandle,'>',$outDir."/merged.tsv") or die "cannot write ".$outDir."/merged.tsv";
	print $mergedHandle MultiInterToStringWithOffsets($dataMultiInter,$chrLenDat);
	close $mergedHandle;
	
	#print string/table $chrLenDat
	open(my $chrLenHandle,'>',$outDir."/merged.chrlen.tsv") or die "cannot write ".$outDir."/merged.chrlen.tsv";
	print $chrLenHandle ChomLengthToString($chrLenDat);
	close $chrLenHandle;
	
	PlotSegData($dataMultiInter,$outDir."/merged.tsv");
	
	warn localtime(time())." [INFO] Done\n";
}

sub GetSourceDir {
	my $srcdir = $0;
	$srcdir =~ s!\w*\.pl$!!g;
	return $srcdir;
}

sub PlotSegData {
	my $dataMultiInter = shift(@_);
	my $mergedData = shift(@_);
	my $sampleCombi;
	for my $sample1 (sort(keys(%{$dataMultiInter->{'samples'}}))){
		for my $sample2 (sort(keys(%{$dataMultiInter->{'samples'}}))){
			
			next if($sampleCombi->{$sample1}->{$sample2} || $sample2 eq $sample1);
			$sampleCombi->{$sample1}->{$sample2}++;
			my $src = GetSourceDir();
			my $cmd = "Rscript $src/PlotSegsSampleVsSample.R $mergedData $sample1 $sample2 2>/dev/null";
			warn "## ".localtime(time())." [INFO] running command for plotting '$cmd'\n";
			warn CmdRunner($cmd);
		}
	}
}

sub ReadDict {
	open( my $DictHandle, '<', $_[0] )
	  or die "cannot read dictionary file '$_[0]' ";
	my %ChrLenData;
	my $RefIndex = 0;
	my $OffSet   = 1;
	while ( my $line = <$DictHandle> ) {
		if ( $line =~ m/^\@SQ/ ) {
			my @tabdelim = split( "\t", $line );
			my @SN       = split( ':',  $tabdelim[1] );
			my @LN       = split( ':',  $tabdelim[2] );
			$ChrLenData{'ChrLen'}{ $SN[1] } = $LN[1];

			$ChrLenData{'ChrOff'}{ $SN[1] } = $OffSet;
			$OffSet = $OffSet + $LN[1];

			$ChrLenData{'order'}{$RefIndex} = $SN[1];
			$RefIndex++;
		}

	}
	$ChrLenData{'len'}=$OffSet;
	$ChrLenData{'REFcount'} = $RefIndex;
	close($DictHandle);
	return \%ChrLenData;
}
sub ChomLengthToString {
	my $self = shift @_;
	my $string = "";
	$string .= "chrom\toffset\tlength\n";
	for my $chrom (sort(keys(%{$self->{'ChrOff'}}))){		
		$string .=$chrom."\t".$self->{'ChrOff'}->{$chrom}."\t".$self->{'ChrLen'}->{$chrom}."\n";
	}
	return $string;
}

sub MultiInter {
	my $bedFiles = shift @_;
	my $dataMultiInter;
	
	for my $bed (@{$bedFiles}){
		my $cmd = "bash -c 'set -e -o pipefail && bedtools multiinter -i ".join(" ", @{$bedFiles})."|cut -f1,2,3 |bedtools intersect -wao -a - -b ".$bed."'";
	
		warn localtime(time())." [INFO] running bedtools for bed annotation as '$cmd'\n";
		open(my $multannothandle,'-|',$cmd);
		while(<$multannothandle>){
			$dataMultiInter=ParseMultiInter($_,$dataMultiInter);
			#warn;
			#die if $. > 10;
		}
		
	}
	#print Dumper($dataMultiInter);
	return $dataMultiInter;
}


sub ParseMultiInter{
	my $l = shift @_;
	my $dataMultiInter=shift @_;
	
	chomp($l);
	my @t=();
	@t = split("\t",$l);
	#my $chrom=$t[0];
	
	my $sample;
	#null intersect== no original annotation in file
	if(($t[6] ne ".") && ($t[6] ne "")){
		%{$sample} = split/;|=/,($t[6]);
		$sample->{'chromSegment'}=$t[3];
		$sample->{'startSegment'}=($t[4]+1);
		$sample->{'endSegment'}=$t[5];
		
		$dataMultiInter->{'samples'}{ $sample->{'ID'} }++;
		
		map{$dataMultiInter->{'annotations'}->{ $_ }++;}(keys(%{$sample}));			
	#												chrom,start0,end
		push(@{$dataMultiInter->{'regions'}->{join("\t",($t[0],($t[1] + 1),$t[2]))}},$sample);
	}
	return $dataMultiInter;
}

sub MultiInterToString{
	my $self=shift @_;
	my $string="";
	
	#header	
	$string.="chrom\tstart\tend";
	for my $sample (sort(keys(%{$self->{'samples'}}))){
	#die Dumper($self->{'samples'});
		for my $annotation (sort(keys(%{$self->{'annotations'}}))){
			$string.="\t".${sample}.'.'.${annotation};
		}
	}
	$string.="\n";
	#die $string;
	#for regions
	for my $region (sort(keys(%{$self->{'regions'}}))){
		#die Dumper($self->{'regions'}->{$region});
		#for samples
		$string.= $region;
		for my $sample (sort(keys(%{$self->{'samples'}}))){
			my $ref='';
			for my $sampledat (@{$self->{'regions'}->{$region}}){
				#die Dumper($self->{'samples'});
				if(defined($sampledat->{'ID'}) && $sampledat->{'ID'} eq $sample){
					$ref=$sampledat;
				}
			}
			for my $annotation (sort(keys(%{$self->{'annotations'}}))){
				 if($ref && defined($ref->{${annotation}})){
				 	$string.="\t".$ref->{${annotation}};
				 }elsif($ref && not(defined($ref->{${annotation}}))){
				 	$string.="\t"; die "Undef error: time to troubleshoot"
				 }else{
				 	$string.="\t";
				 }
			}
		}
		$string.="\n";
	}
	print "$string";
	die "$string";
	
}

sub MultiInterToStringWithOffsets{
	my $self=shift @_;
	my $offsets=shift @_;
	
	my $string="";
	
	#header	
	$string.="chrom\tstart\tend\tstart.offset\tend.offset";
	for my $sample (sort(keys(%{$self->{'samples'}}))){
	#die Dumper($self->{'samples'});
		for my $annotation (sort(keys(%{$self->{'annotations'}}))){
			$string.="\t".${sample}.'.'.${annotation};
		}
	}
	$string.="\n";
	#die $string;
	#for regions
	for my $region (sort(keys(%{$self->{'regions'}}))){
		#die Dumper($self->{'regions'}->{$region});
		#for samples
		$string.= $region;
		my ($chrom,$start, $end) = split('\t',$region);
		my $offset=$offsets->{'ChrOff'}->{$chrom};
		$string.="\t".($offset+$start)."\t".($offset+$end);
		
		for my $sample (sort(keys(%{$self->{'samples'}}))){
			my $ref='';
			for my $sampledat (@{$self->{'regions'}->{$region}}){
				#die Dumper($self->{'samples'});
				if(defined($sampledat->{'ID'}) && $sampledat->{'ID'} eq $sample){
					$ref=$sampledat;
				}
			}
			for my $annotation (sort(keys(%{$self->{'annotations'}}))){
				 if($ref && defined($ref->{${annotation}})){
				 	$string.="\t".$ref->{${annotation}};
				 }elsif($ref && not(defined($ref->{${annotation}}))){
				 	$string.="\t"; die "Undef error: time to troubleshoot"
				 }else{
				 	$string.="\t";
				 }
			}
		}
		$string.="\n";
	}
	#print "$string";
	#die "$string";
	return  "$string";
}

sub SegNew{
	my %init = @_;
	my $self;
	if(-e $init{'file'}){
		
		$self=SegOpen(\%init);
	}else{
		die "SegNew:needs segfile! as : SegNew(file=>'file.seg')";
	}
	return $self;
}
sub SegOpen {
	my $self = shift @_;
	#print Dumper($self);
	if($$self{'file'} && not($$self{'handle'})){
		open($$self{'handle'},'<',$$self{'file'}) or die "SegOpen:cannot read segfile!specify as : SegNew(file=>'file.seg')";
	}else{
		die "SegOpen:needs segfile! as : SegNew(file=>'file.seg')";
	}
	return $self;
}

sub SegClose {
	my $self=shift @_;
	close $$self{'handle'} or die 'fail on close of segfile'.$$self{'file'};
	#                          ^ warn or die?
	return;
}

sub SegNextDataHash{
	my $self = shift @_;
	if($self->{'handle'}){
		return if(eof($self->{'handle'}));
		
		$self->{'record'}->{'line'} = readline($self->{'handle'});
		
		if($. == 1){
			#parseheader here ?
			$self = SegParseHeader($self);
			$self->{'record'}->{'line'} = readline($self->{'handle'});	
		}
		$self=SegParseLine($self);
	}else{
		die "SegNextDataHash:needs seghandle!".$$self{'handle'};
	}
	return $self;
}

sub SegParseLine{
	my $self = shift @_;
	if($self->{'record'}->{'line'}){
		chomp($self->{'record'}->{'line'});
		@{$self->{'record'}->{'array'}}=split("\t",$self->{'record'}->{'line'});
		
		
		$self->{'record'}->{'chrom'}= $self->{'record'}->{'array'}->[1];
		$self->{'record'}->{'start'}= $self->{'record'}->{'array'}->[2];
		$self->{'record'}->{'end'}= $self->{'record'}->{'array'}->[3];
		
		my %sample=(	'ID'		=> $self->{'record'}->{'array'}->[0],
						'num.mark'	=> $self->{'record'}->{'array'}->[4],
						'seg.mean'	=> $self->{'record'}->{'array'}->[5]);		
		@{$self->{'record'}->{'samples'}}=();
		push(@{$self->{'record'}->{'samples'}},\%sample);
		
	}else{
		die "SegParseLine:needs segline!".$self->{'line'};
	}
	return $self;
}
sub SegParseHeader {
	my $self = shift @_;
	if($$self{'record'}{'line'}){
		#do something
		chomp($self->{'record'}->{'line'});
		warn localtime(time())." [INFO] Skipping first line / header parse: '".$self->{'record'}->{'line'}."'";
	}else{
		die "SegParseLine:needs segline!".$self->{'line'};
	}
	return $self;	
}
sub SegDataHashAsBedString{
	my $self = shift @_;
	my $bedString = join("\t",(	$$self{'record'}{'chrom'},
								($$self{'record'}{'start'} - 1),
								$$self{'record'}{'end'},
								'ID='. $self->{'record'}->{'samples'}->[0]->{'ID'}.';num.mark='. $self->{'record'}->{'samples'}->[0]->{'num.mark'}.';seg.mean='. $self->{'record'}->{'samples'}->[0]->{'seg.mean'}))."\n";
	return $bedString;
}

sub CmdRunner {
	my $ret;
	my $cmd = join(" ",@_);
	
	warn localtime( time() ). " [INFO] system call:'". $cmd."'.\n";
	
	@{$ret} = `($cmd )2>&1`;
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

	return @{$ret};
}
