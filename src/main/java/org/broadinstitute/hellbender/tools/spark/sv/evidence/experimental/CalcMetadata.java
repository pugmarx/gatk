package org.broadinstitute.hellbender.tools.spark.sv.evidence.experimental;

import org.apache.logging.log4j.LogManager;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.SVReadFilter;

import java.util.Collections;

@CommandLineProgramProperties(summary="Dump some data about the reads.",
        oneLineSummary="Dump some data about the reads.",
        programGroup=StructuralVariationSparkProgramGroup.class)
@BetaFeature
public class CalcMetadata extends GATKSparkTool {
    private static final long serialVersionUID = 1L;


    @Override public boolean requiresReads() { return true; }

    @Override protected void runTool( final JavaSparkContext ctx ) {
        final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params =
                new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();
        final ReadMetadata readMetadata =
                new ReadMetadata(Collections.emptySet(), getHeaderForReads(),
                        10000, getUnfilteredReads(),
                        new SVReadFilter(params), LogManager.getLogger(CalcMetadata.class));
        System.out.println("Effective Kmer coverage for K=" + params.kSize + " = " +
                readMetadata.getAccurateKmerCoverage(params.kSize));
        ReadMetadata.writeMetadata(readMetadata, "/dev/stdout");
    }
}
