use 5.020002;
use ExtUtils::MakeMaker;


my $extdeps;
$extdeps = [
	{
		"dep" => "R",
		"test" => "R --version",
		"hint" => "wget https://cran.r-project.org/src/base/R-3/R-3.4.4.tar.gz && tar -zxf R* && cd R* && configure && make" 
		#maybe install?
	},{	"dep" => "htslib",
		"test" => "(exit 0)",
		"hint" => "wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && tar -jxf htslib* && cd htslib* && make" 
	},{	"dep" => "samtools",
		"test" => "which samtools",
		"hint" => "wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && tar -jxf samtools* && cd samtools* && make"
	},{	"dep" => "bedtools",
		"test" => "bedtools --version",
		"hint" => "wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && tar -xzf https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && cd bed* && make"
	},{	"dep" => "bwa",
		"test" => "which bwa", # the actual command exits with 1 with no working args
		"hint" => "wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && tar -xjf https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && cd bwa-* && make "
	},{	"dep" => "DNAcopy the R library",
		"test" => "Rscript -e 'library(\"DNAcopy\")'",
		"hint" => "Rscript -e 'install.packages(c(\"DNAcopy\"))'"
	},{	"dep" => "tabix",
		"test" => "which tabix",
		"hint" => "sudo apt get install tabix"
        },{	"dep" => "vcftools",
		"test" => "vftools",
		"hint" => "sudo apt get install vcftools"
	}
];

for my $dep (@{$extdeps}){
	warn "testing for ".$dep -> {'dep'}.".\n";
	CmdRunner($dep -> {'test'});
}


sub CmdRunner {
        my $ret;
        my $cmd = join(" ",@_);

        warn localtime( time() ). " [INFO] system call:'". $cmd."'.\n";

        @{$ret} = `($cmd )2>&1`;
	my $exitcode = $?;
        if ($exitcode == -1) {
                die localtime( time() ). " [ERROR] failed to execute: $!\n";
        }
	elsif ($exitcode & 127) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($exitcode & 127),  ($exitcode & 128) ? 'with' : 'without';
        }
	elsif ($exitcode != 0) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($exitcode & 127),  ($exitcode & 128) ? 'with' : 'without';
        }else {
               	warn localtime( time() ). " [INFO] " . sprintf "child exited with value %d\n", $exitcode >> 8;
        }
		return @{$ret};
}

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
	"NAME"              => 'pipeline::util',
	"VERSION_FROM"      => 'lib/pipeline/util.pm', # finds $VERSION, requires EU::MM from perl >= 5.5
	"PREREQ_PM"         => {
		#'Vcf' => '0.1.12b', #from vcftools
		'List::MoreUtils' => '0.416',
	}, # e.g., Module::Name => 1.1
	"BUILD_REQUIRES" => {
		'Test::More' => '0.47',
		'ExtUtils::MakeMaker' => '0'
	},
	"EXE_FILES" => [
		'src/AddAdBasedAnnotations.pl',
		'src/AddAlleleFrequenciesToSeg.pl',
		'src/CalleriseVcf.pl',
		'src/CollectNugeneLandingProbeMetrics.pl',
		'src/dbNSFPAnnotator.pl',
		'src/DumpFieldsForVariantsToTable.pl',
		'src/filterCombinedVariantsForGatk.pl',
		'src/FilterFreeBayes.pl',
		'src/ForGerardTeMeermanVcfAnnotation.pl',
		'src/GenerateTableDescriptionByVcfHeader.pl',
		'src/GTfGet1000bpExonsBeforeTES.pl',
		'src/multiIntersectSeg.pl',
		'src/MutectAnnotationsToSampleFormat.pl',
		'src/NugeneDigitalSplitter.pl',
		'src/PlotFloatsOnInterVals.pl',
		'src/RecoverSampleAnnotationsAfterCombineVariants.pl',
		'src/refineSam.pl',
		'src/SampleSheetAddVal.pl',
		'src/tickerRefine.pl',
		'src/tickertape.pl',
		'src/trimByBed.pl',
		'src/updateGTPLByAD.pl',
		'src/updateSimilarPos.pl',
		'src/VcfaddAdBasedZScores.pl',
		'src/VcfSnpEffAsGatk.pl',
	],

	"ABSTRACT_FROM"     => 'lib/pipeline/util.pm', # retrieve abstract from module
	"AUTHOR"            => 'm.m.terpstra',
	#LICENSE           => 'perl',
	#Value must be from legacy list of licenses here
	#http://search.cpan.org/perldoc?Module%3A%3ABuild%3A%3AAPI
);
