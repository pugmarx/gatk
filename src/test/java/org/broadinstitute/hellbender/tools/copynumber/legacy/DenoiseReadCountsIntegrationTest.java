package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Integration tests for {@link DenoiseReadCounts}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class DenoiseReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/legacy/coverage/denoising";
    private static final File WGS_READ_COUNTS_TSV_FILE = new File(TEST_SUB_DIR, "denoise-read-counts-wgs-read-counts-HCC1143_BL-n1-chr20-downsampled-deduplicated.readCounts.tsv");
    private static final File WGS_READ_COUNTS_HDF5_FILE = new File(TEST_SUB_DIR, "denoise-read-counts-wgs-read-counts-HCC1143_BL-n1-chr20-downsampled-deduplicated.readCounts.hdf5");
    private static final File WGS_ANNOTATED_INTERVALS_FILE = new File(TEST_SUB_DIR, "denoise-read-counts-wgs-annotated-intervals.annotated.tsv");
    private static final File WGS_NO_GC_PON_FILE = new File(largeFileTestDir, "cnv_somatic_workflows_test_files/wgs-no-gc.pon.hdf5");
    private static final File WGS_DO_GC_PON_FILE = new File(largeFileTestDir, "cnv_somatic_workflows_test_files/wgs-do-gc.pon.hdf5");
    private static final List<Integer> NUMBER_OF_EIGENVALUES_LIST = Arrays.asList(null, 1, 10);

    //create all combinations of arguments
    @DataProvider(name = "dataDenoiseReadCounts")
    public Object[][] dataDenoiseReadCounts() {
        final List<List<Object>> data = new ArrayList<>();
        for (final File inputReadCountsFile : Arrays.asList(WGS_READ_COUNTS_TSV_FILE, WGS_READ_COUNTS_HDF5_FILE)) {
            for (final File annotatedIntervalsFile : Arrays.asList(WGS_ANNOTATED_INTERVALS_FILE, null)) {
                for (final File ponFile : Arrays.asList(WGS_NO_GC_PON_FILE, WGS_DO_GC_PON_FILE, null)) {
                    for (final Integer numberOfEigenvalues : NUMBER_OF_EIGENVALUES_LIST) {
                        final ArgumentsBuilder arguments = new ArgumentsBuilder()
                                .addFileArgument(StandardArgumentDefinitions.INPUT_SHORT_NAME, inputReadCountsFile);
                        if (annotatedIntervalsFile != null) {
                            arguments.addFileArgument(CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_SHORT_NAME, annotatedIntervalsFile);
                        }
                        if (ponFile != null) {
                            arguments.addFileArgument(CopyNumberStandardArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, ponFile);
                        }
                        if (numberOfEigenvalues != null) {
                            arguments.addArgument(CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_SHORT_NAME, numberOfEigenvalues.toString());
                        }
                        data.add(Arrays.asList(arguments, ponFile == null));    //set isStandardizedEqualsDenoised = true if no PoN
                    }
                }
            }
        }
        return data.stream().map(List::toArray).toArray(Object[][]::new);
    }

    @Test(dataProvider = "dataDenoiseReadCounts")
    public void testDenoiseReadCounts(final ArgumentsBuilder argumentsBuilder,
                                      final boolean isStandardizedEqualsDenoised) {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = argumentsBuilder
                .addFileArgument(CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile)
                .addFileArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile)
                .addArgument(StandardArgumentDefinitions.VERBOSITY_NAME, "INFO")
                .getArgsArray();
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names and intervals should always be the same
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        //standardized and denoised copy ratios should be the same if PoN is not provided
        Assert.assertEquals(standardizedCopyRatios.getLog2CopyRatioValues().equals(denoisedCopyRatios.getLog2CopyRatioValues()), isStandardizedEqualsDenoised);
    }
}