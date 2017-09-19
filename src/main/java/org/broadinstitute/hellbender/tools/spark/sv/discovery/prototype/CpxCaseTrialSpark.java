package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.FileUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.SvDiscoverFromLocalAssemblyContigAlignmentsSpark.writeSAM;

@CommandLineProgramProperties(summary="Parses a SAM file containing long reads alignments, and outputs cxSV rearrangements.",
        oneLineSummary="Parses a long read SAM file, and outputs cxSV rearrangements.",
        usageExample = "CpxCaseTrialSpark \\" +
                "-I ~/Desktop/summary/Cpx.sam \\" +
                "-O ~/Desktop/cpxCases",
        omitFromCommandLine = true,
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class CpxCaseTrialSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger();

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference, human reference (hg19 or hg38) assumed when omitted",
            shortName = "nonCanoChrFile",
            fullName = "nonCanonicalChromosomeNamesFile", optional = true)
    private String nonCanonicalChromosomeNamesFile;

    @Argument(doc = "output directory", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getUnfilteredReads();
        final SAMFileHeader header = getHeaderForReads();
        final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(header);

        final JavaRDD<AlignedContig> contigsWithAlignmentsReconstructed =
                FilterLongReadAlignmentsSAMSpark.filterByScore(reads, header, nonCanonicalChromosomeNamesFile, localLogger)
                        .filter(lr -> lr.alignmentIntervals.size()>1).cache();

        final JavaRDD<AlignedContig> cpxContigs =
                SvDiscoverFromLocalAssemblyContigAlignmentsSpark
                        .divertReadsByPossiblyRawTypes(contigsWithAlignmentsReconstructed, localLogger)
                        .get(SvDiscoverFromLocalAssemblyContigAlignmentsSpark.RawTypes.Cpx).cache();

        localLogger.info(cpxContigs.count() + " contigs suggesting complex variants.");

        FileUtils.createDirToWriteTo(outputDir);
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> split = RDDUtils.split(cpxContigs, tig -> tig.alignmentIntervals.size() == 2, false);

        // split by ==2 alignments
        final JavaRDD<AlignedContig> twoMappings = split._1;
        final JavaRDD<AlignedContig> moreThanTwoMappings = split._2.cache();

        localLogger.info(twoMappings.count() + " contigs have 2 alignments after filtering.");
        writeSAM(twoMappings, "twoMappings", reads, headerBroadcast, outputDir, null);
        cpxContigs.unpersist();

        // split by if HeadChr==TailChr
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> split1 = RDDUtils.split(moreThanTwoMappings,
                tig -> {
                    final AlignmentInterval head = tig.alignmentIntervals.get(0);
                    final AlignmentInterval tail = tig.alignmentIntervals.get(tig.alignmentIntervals.size() - 1);
                    return head.referenceSpan.getContig().equals(tail.referenceSpan.getContig());
                },
                true);
        final JavaRDD<AlignedContig> cameBack = split1._1.cache();
        final JavaRDD<AlignedContig> didnotCameBack = split1._2.cache();
        localLogger.info(didnotCameBack.count() + " contigs have more than 2 mappings that have head and tail mapped to diff chr.");
        localLogger.info(cameBack.count() + " contigs have more than 2 mappings that have head and tail mapped to the same chr.");
        writeSAM(didnotCameBack, "2+DidNotCameBack", reads, headerBroadcast, outputDir, null);

        // split by if HeadStrand==TailStrand given HeadChr==TailChr
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> split2 = RDDUtils.split(cameBack,
                tig -> {
                    final AlignmentInterval head = tig.alignmentIntervals.get(0);
                    final AlignmentInterval tail = tig.alignmentIntervals.get(tig.alignmentIntervals.size() - 1);
                    return head.forwardStrand == tail.forwardStrand;
                }, false);
        final JavaRDD<AlignedContig> sameStrand = split2._1;
        final JavaRDD<AlignedContig> diffStrand = split2._2;
        localLogger.info(sameStrand.count() + " contigs have head/tail mapped to the same chr with the same strand");
        localLogger.info(diffStrand.count() + " contigs have head/tail mapped to the same chr with opposite strands");
        writeSAM(sameStrand, "2+SameChrSameStrand", reads, headerBroadcast, outputDir, null);
        writeSAM(diffStrand, "2+SameChrDiffStrand", reads, headerBroadcast, outputDir, null);
    }
}
