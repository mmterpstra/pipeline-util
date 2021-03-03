use 5.020002;
use Data::Dumper;
use ExtUtils::MakeMaker;

#How to add a new file
# -create the script and put in into src
# -update the tests to do al least a syntax test
# -Update pipeline/util with newer version numer
#git add & commit it


my $extdeps;
$extdeps = [
	{
		"dep" => "R",
		"test" => "command -v -- R",
		"hint" => "wget https://cran.r-project.org/src/base/R-3/R-3.4.4.tar.gz && tar -zxf R* && cd R* && configure && make" 
		#maybe install?
	},{	"dep" => "htslib",
		#"test" => "command -v -- htsfile",
		#"test"=> "pkg-config --libs  htslib || pkg-config --libs  htslib1",
		"test" => "exit 0",
		"hint" => "wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && tar -jxf htslib* && cd htslib* && make" 
	},{	"dep" => "samtools",
		"test" => "command -v -- samtools",
		"hint" => "wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && tar -jxf samtools* && cd samtools* && make"
	},{	"dep" => "bedtools",
		"test" => "command -v -- bedtools",
		"hint" => "wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && tar -xzf https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && cd bed* && make"
	},{	"dep" => "bwa",
		"test" => "command -v -- bwa", # the actual command exits with 1 with no working args
		"hint" => "wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && tar -xjf https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && cd bwa-* && make "
	},{	"dep" => "DNAcopy the R library",
		"test" => "Rscript -e 'library(\"DNAcopy\")'",
		"hint" => "Rscript -e 'install.packages(c(\"DNAcopy\"))'"
	},{	"dep" => "tabix",
		"test" => "command -v -- tabix",
		"hint" => "sudo apt get install tabix"
        },{	"dep" => "vcftools",
		"test" => "command -v -- vcftools",
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
	
	#adding a `set -e; set -o pipefail` although safer wil crash the travis build env. 
        @{$ret} = `( $cmd )2>&1`;
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

my @exe_files = glob("src/*.pl");
push(@exe_files,glob("src/*.Rscript"));

my $pm = 'lib/pipeline/util.pm';

my $git;
@{$git} = CmdRunner("which git || exit 0");
warn Dumper(\$git);
if($git -> [0] =~ m!.*/git!){
	#git present
	warn localtime( time() ). " [INFO] git present";
	my $git_descr;
	@{$git_descr} = CmdRunner("git describe --tag || exit 0");
	if($git_descr -> [0] =~ m/\d+\.\d+\.\d+/){ 
		warn localtime( time() ). " [INFO] in git repo with tag ".join(" ",@{$git_descr}); 
		#this likely means were in a repo so we can edit the $VERSION if these arent the same as the github repo
		my $ver_tag; @{$ver_tag} = CmdRunner("git describe --tags| perl -wpe 's/\n//'");
		my $ver_pm;@{$ver_pm} = CmdRunner("grep '\$VERSION' $pm | perl -wpe 's/.*=\\s*\"|\";\\n//g'");
		warn localtime( time() ). " [INFO] " . " tags detectected git: " . $ver_tag -> [0] . " pm tag: ".$ver_pm -> [0]."\n";
		if($ver_tag -> [0] ne $ver_pm -> [0] ){
			warn  localtime( time() ). " [INFO] Updating version to match with git tag";
			my $new_tag = $ver_tag -> [0];
			my $backup_tag = $ver_pm -> [0];
			CmdRunner("perl -i${backup_tag}.bak -wpe 's/\^our \\\$VERSION.*/our \\\$VERSION = \"$new_tag\";/' $pm");
		}
	}

}

WriteMakefile(
	"NAME"              => 'pipeline::util',
	"VERSION_FROM"      => $pm, # finds $VERSION, requires EU::MM from perl >= 5.5
	"PREREQ_PM"         => {
		#'Vcf' => '0.1.12b', #from vcftools
		'List::MoreUtils' => '0.416',
	}, # e.g., Module::Name => 1.1
	"BUILD_REQUIRES" => {
		'Test::More' => '0.47',
		'ExtUtils::MakeMaker' => '0'
	},
	"EXE_FILES" => [@exe_files],
	"ABSTRACT_FROM"     => 'lib/pipeline/util.pm', # retrieve abstract from module
	"AUTHOR"            => 'm.m.terpstra',
	#LICENSE           => 'perl',
	#Value must be from legacy list of licenses here
	#http://search.cpan.org/perldoc?Module%3A%3ABuild%3A%3AAPI
);
