# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl pipeline-util.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;


my $exes = [
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
	'src/VcfQssFix.pl',
	'src/VcfTableExportOneVariantPerSample.pl',
];


use Test::More tests => 30;
BEGIN { use_ok('pipeline::util') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

for my $exe (@{$exes}){
	ok(system('bash -c "echo > /dev/stderr && set -ex && perl -wc '.$exe.'"') == 0, "$exe syntax test");
}
# test bedtrim
ok(system( 'bash -c "echo > /dev/stderr && set -e && ' .
	'perl src/trimByBed.pl -s t/data/bedtrim_before.sam  -b t/data/bedtrim_probe.bed -o t/data/bedtrim_test -n t/data/bedtrim_test && ' . 
	'diff <(samtools view t/data/bedtrim_after.bam ) <(samtools view t/data/bedtrim_test.bam )"' ) == 0, "bedtrim minimal functional test");

ok(system('bash -c "echo > /dev/stderr && set -e -o pipefail && '.
	'diff <(export PERL5LIB=\"blib/lib/\":$PERL5LIB && perl src/RecoverSampleAnnotationsAfterCombineVariantsByPosWalk.pl complex.vcf t/data/annot.vcf t/data/annot.call1.vcf t/data/annot.call2.vcf 2>/dev/null) '.
	' t/data/recov.vcf &>/dev/stderr"') == 0 , 'RecoverSampleAnnotationsAfterCombineVariantsByPosWalk funtional test');
