package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class PlotDenoisedCopyRatiosIntegrationTest extends CommandLineProgramTest {
    private static String TEST_SUB_DIR = toolsTestDir + "copynumber/legacy/plotting/";

    //test files
    private static final File STANDARDIZED_COPY_RATIOS_FILE = new File(TEST_SUB_DIR, "plotting-copy-ratios.tsv"); //just use the TN file
    private static final File DENOISED_COPY_RATIOS_FILE = new File(TEST_SUB_DIR, "plotting-copy-ratios.tsv");
    private static final File SEQUENCE_DICTIONARY_FILE = new File(TEST_SUB_DIR, "plotting-sequence-dictionary.dict");

    //test files for invalid configurations
    private static final File SEQUENCE_DICTIONARY_WITH_NO_CONTIGS_ABOVE_MINIMUM_LENGTH_FILE = new File(TEST_SUB_DIR, "plotting-sequence-dictionary-with-no-contigs-above-minimum-length.dict");
    private static final File COPY_RATIOS_WITH_BAD_SAMPLE_NAME_FILE = new File(TEST_SUB_DIR, "plotting-copy-ratios-with-bad-sample-name.tsv");
    private static final File COPY_RATIOS_OUT_OF_BOUNDS_FILE = new File(TEST_SUB_DIR, "plotting-copy-ratios-out-of-bounds.tsv");

    private static final String OUTPUT_PREFIX = "test";
    private static final int THRESHOLD_PLOT_FILE_SIZE_IN_BYTES = 50000;  //test that data points are plotted (not just background/axes)

    //checks that output files with reasonable file sizes are generated, but correctness of output is not checked
    @Test()
    public void testPlotting() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, STANDARDIZED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_Before_After_CR_Lim_4.png").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_Before_After_CR_Lim_4.png").length() > THRESHOLD_PLOT_FILE_SIZE_IN_BYTES);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_Before_After.png").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_Before_After.png").length() > THRESHOLD_PLOT_FILE_SIZE_IN_BYTES);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_preQc.txt").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_preQc.txt").length() > 0);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_postQc.txt").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_postQc.txt").length() > 0);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_dQc.txt").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_dQc.txt").length() > 0);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_scaled_dQc.txt").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_scaled_dQc.txt").length() > 0);
        final double postQc = ParamUtils.readValuesFromFile(new File(outputDir, OUTPUT_PREFIX + "_preQc.txt"))[0];
        final double preQc = ParamUtils.readValuesFromFile(new File(outputDir, OUTPUT_PREFIX + "_postQc.txt"))[0];
        final double dQc = ParamUtils.readValuesFromFile(new File(outputDir, OUTPUT_PREFIX + "_dQc.txt"))[0];
        final double scaled_dQc = ParamUtils.readValuesFromFile(new File(outputDir, OUTPUT_PREFIX + "_scaled_dQc.txt"))[0];
        Assert.assertEquals(dQc, preQc - postQc);
        Assert.assertEquals(scaled_dQc, (preQc - postQc) / preQc);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMinimumContigLength() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, STANDARDIZED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_WITH_NO_CONTIGS_ABOVE_MINIMUM_LENGTH_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testOutputDirExists() throws IOException {
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, STANDARDIZED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,  DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "Non-existent-path",
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingStandardizedFile() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, "Non-existent-file",
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,  DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingDenoisedFile() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, STANDARDIZED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,  "Non-existent-file",
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingSequenceDictionaryFile() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, STANDARDIZED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME,  "Non-existent-file",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testStandardizedSampleNameMismatch() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, COPY_RATIOS_WITH_BAD_SAMPLE_NAME_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testStandardizedDataOutOfBounds() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, COPY_RATIOS_OUT_OF_BOUNDS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDenoisedSampleNameMismatch() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, STANDARDIZED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, COPY_RATIOS_WITH_BAD_SAMPLE_NAME_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDenoisedDataOutOfBounds() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, STANDARDIZED_COPY_RATIOS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, COPY_RATIOS_OUT_OF_BOUNDS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }
}