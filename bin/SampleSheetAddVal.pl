use warnings;
use strict;
use Data::Dumper;
my $use =<<"END";
use: perl $0 samplesheet.csv  prefix> samplesheet.filenames.csv 
END
#@ARGV=("/media/terpstramm/My_Book/rawdata/0549_FFPE_NuGenOnco/Samplesheet_run0549.csv","/media/terpstramm/My_Book/rawdata/0549_FFPE_NuGenOnco/","/media/terpstramm/My_Book/rawdata/tmp01/");

main();
#help
#internalId,project,sampleName,samplePrep,sequencer,seqType,sequencerId,flowcellId,run,lane,barcode,reads1FqGz,reads2FqGz
#internalSampleID,lane,sequencer,Sample,sequencingStartDate,run,flowcell,seqType,barcodeType,barcode,externalSampleID,"Merge with",project,contact,"Sample Type",arrayFile,arrayID,capturingKit,prepKit,GAF_QC_Name,GAF_QC_Date,GAF_QC_Status,GCC_Analysis,"Plates in stock - DNA sequencing (whole genome / capturing)","Rejected for processing","Plates in stock - RNA sequencing","Barcode 2 "

sub main {
	my $samplesheetfile = shift @ARGV;
	my $addition = "=,=";
	$addition = shift @ARGV;#should be "sample=foo,project=bar,key=val"
	warn "## ".localtime(time())." ## INFO ## $0 init with samplesheetfile='$samplesheetfile';addition='$addition'.\n";
	$addition = ParseAddition($addition);
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	#die Dumper($samplesheet);
	$samplesheet = AddVal($samplesheet,$addition);
	#die Dumper($samplesheet);
	print SamplesheetAsString($samplesheet);
}

sub ParseAddition {
	my $string = shift @_;
	my $hash;
	die "prefix should be formatted as 'project=projectname,sampleName=sample1,count=1000000'." if(not($string =~ m/\S*\=\S*\,\S*\=\S*\,\S*\=\S*/));
	for my $pair (CommaseparatedSplit($string)){
		my ($key,$val) = split('=',$pair);
		$hash -> {$key} = $val;
	}
	SelfRequire(%{$hash}, 'req'=> ['sampleName', 'project']);
	return $hash;
}

sub ReadSamplesheet {
	my $samplesheetf = shift @_;
	my $samplesheet;
	open( my $samplesheeth,'<', $samplesheetf) or die "Cannot read samplesheet file:'$samplesheetf'";
	$_=<$samplesheeth>;
	chomp;
	my @h = CommaseparatedSplit($_);
	#die Dumper(\@h);
	while(<$samplesheeth>){
		chomp;
		my @c = CommaseparatedSplit($_);
		#die Dumper(\@c);
		if(scalar(@c)==scalar(@h)){
			my %d;
			my $i=0;
			map{$d{$_}=$c[$i]; $i++}(@h);
			$c[$i]=join(",",@h);
			#ReadFileNameConstructor(\%d);
			push(@$samplesheet,\%d);
		}else{
			die "Header is not of equal length compared sample line";
		}
		#$d{'samplesheetSampleNo'}=
	}
	return $samplesheet;
}
sub CommaseparatedSplit {
	my $string=pop @_;
	#needs to be fixed for citation marks!
	warn "Line contains citation marks: this is currently not supported!!. I hope this works. Line=$_" if($string =~ /"|'/);
	my $i = index($string,",");
	if( $i > -1){
		push(@_,substr($string,0,$i));
		push(@_,substr($string,$i+1));
		@_ = CommaseparatedSplit( @_ );
	}else{
		push(@_,$string);
		return @_;
	}
}

sub SamplesheetAsString {
	my $self = shift @_;
	my $string = '';
	#get header values;
	my %h;
	for my $sample (@$self){
		for my $key (keys(%$sample)){
			$h{$key}++;
		}
	}
	my @h = sort {$b cmp $a} (keys(%h));
	$string.=join(",",@h)."\n";
	#warn scalar(@$self);
	for my $sample (@$self){
		my @c;
		for my $h (@h){
			push (@c,$$sample{$h});
		}
		$string.=join(",",@c)."\n";
	}
	return $string;
}

sub SelfRequire {
        #calls like SelfRequire(%{$self}, 'req'=> ['vcf', 'fields']);
        my $self;
        %{$self}= @_;
        #die Dumper($self)." ";
        my $reqs = $self -> {'req'} or die "no requirement array as input";
        for my $tag (@{$reqs}){
                defined($self -> {$tag}) or die "no requirement '$tag' as input";
        }
}


sub AddVal {
	my $samplesheet = shift @_;
	my $samplesheetNew;
	my $additiondescr = shift @_;
	for my $sample (@{$samplesheet}){
		if ($additiondescr -> {'project'} eq $sample -> {'project'} && $additiondescr -> {'sampleName'} eq $sample -> {'sampleName'} ){
			for my $key (keys( %{$additiondescr})){
				next if($key eq 'project' ||  $key eq 'sampleName');
				 $sample -> { $key } = $additiondescr -> {$key};
			}
		}
		push(@{$samplesheetNew}, $sample)
	}
	return $samplesheetNew;
}
