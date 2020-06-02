#!/usr/bin/env perl	
use warnings;
use strict;
use Data::Dumper;
use Carp;
main();

sub main{
	warn "Commandline: $0 ".join(' ',@ARGV)."\nUse: $0 chrcol rename.tsv input.tsv > output.tsv\n rename.tsv content:chromsome to convert{tab}new chromosome name{newline}";
	my $chrcol = shift(@ARGV);
	my $renametsv = shift(@ARGV);
	my $inputtsv = shift(@ARGV);
	#process rename.tsv
	my $rename;#example below
	#$VAR2 = {
	#          '1' => {
	#                   'to' => 'chr1',
	#                   'from' => '1',
	#                   'line' => '1'
	#                 },
	#          'GL00220' => {
	#                         'to' => 'chrUn_GL00220',
	#                         'from' => 'GL00220',
	#                         'line' => '3'
	#                       },
	#          '2' => {
	#                   'to' => 'chr2',
	#                   'from' => '2',
	#                   'line' => '2'
	#                 }
	#        };
	open(my $renamehandle,'<',$renametsv) or die "cannot read $renametsv. Is the file ok?";
	while(<$renamehandle>){
		next if ($_ eq "\n");
		my $record;
		$record = DumbReader($_);
		$rename -> {$record -> [0]} = {"from" => $record -> [0],"to" => $record -> [1],"line" => $.};
		#warn Dumper($record,$rename)."";
	}
	warn "## INFO ## Reading input file succesful rename info: " . Dumper($rename);
	close($renamehandle);
	#process bed
	
	open(my $inputhandle,'<',$inputtsv) or die "cannot read $inputtsv. Is the file ok?";
	
	while(<$inputhandle>){
		my $record;
		$record = DumbReader($_);
		#faster than regex
		if(substr($record -> [0],0,1) eq '#'){
			#print and next on comment lines
			print RecordAsString($record);
		}
		elsif(defined($rename -> {$record -> [$chrcol]})){
			#warn "Oldrecord: ".RecordAsString($record);
			$record -> [$chrcol] = $rename -> {$record -> [$chrcol]} -> {'to'};
			print RecordAsString($record);
			
		}else{
			warn "Omitted the following region due to chromosome not present in rename.tsv : '".Dumper($record)."'";
		}
	}
	close($inputhandle);
}

sub DumbReader{
	$_= shift(@_);
	chomp;
	my $record;
	@{$record}=split("\t");
	return $record;
}
sub RecordAsString{
	my $r= shift(@_);
	#warn Dumper($_);
	return join("\t",@{$r})."\n";
}

