package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * SparkTool to identify 63-mers in the reference that occur more than 3 times.
 */
@CommandLineProgramProperties(summary="Find the set of high copy number kmers in a reference.",
        oneLineSummary="find ref kmers with high copy number",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class FindBadGenomicKmersSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    @VisibleForTesting static final int MAX_KMER_FREQ = 3;

    @Argument(doc = "file for ubiquitous kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "kmer size", fullName = "kSize", optional = true)
    private int kSize = StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.KMER_SIZE;

    @Argument(doc = "maximum kmer DUST score", fullName = "kmerMaxDUSTScore")
    private int maxDUSTScore = StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.MAX_DUST_SCORE;

    @Argument(doc = "additional high copy kmers (mitochondrion, e.g.) fasta file name",
            fullName = "highCopyFasta", optional = true)
    private String highCopyFastaFilename;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /** Get the list of high copy number kmers in the reference, and write them to a file. */
    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader hdr = getHeaderForReads();
        SAMSequenceDictionary dict = null;
        if ( hdr != null ) dict = hdr.getSequenceDictionary();
        final ReferenceMultiSource referenceMultiSource = getReference();
        Collection<SVKmer> killList = findBadGenomicKmers(ctx, kSize, maxDUSTScore, referenceMultiSource, dict);
        if ( highCopyFastaFilename != null ) {
            killList = SVUtils.uniquify(killList, processFasta(kSize, maxDUSTScore, highCopyFastaFilename));
        }

        SVUtils.writeKmersFile(kSize, outputFile, killList);
    }

    /** Find high copy number kmers in the reference sequence */
    @VisibleForTesting
    static List<SVKmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                             final int kSize,
                                             final int maxDUSTScore,
                                             final ReferenceMultiSource ref,
                                             final SAMSequenceDictionary readsDict ) {
        // Generate reference sequence RDD.
        final JavaRDD<byte[]> refRDD = SVUtils.getRefRDD(ctx, kSize, ref, readsDict,
                                                        SVReferenceUtils.REF_RECORD_LEN, SVReferenceUtils.REF_RECORDS_PER_PARTITION);

        // Find the high copy number kmers
        return SVReferenceUtils.processRefRDD(kSize, maxDUSTScore, MAX_KMER_FREQ, refRDD);
    }

    @VisibleForTesting static List<SVKmer> processFasta( final int kSize,
                                                         final int maxDUSTScore,
                                                         final String fastaFilename) {
        try ( BufferedReader rdr = new BufferedReader(new InputStreamReader(BucketUtils.openFile(fastaFilename))) ) {
            final List<SVKmer> kmers = new ArrayList<>((int) BucketUtils.fileSize(fastaFilename));
            String line;
            final StringBuilder sb = new StringBuilder();
            final SVKmer kmerSeed = new SVKmerLong();
            while ( (line = rdr.readLine()) != null ) {
                if ( line.charAt(0) != '>' ) sb.append(line);
                else if ( sb.length() > 0 ) {
                    SVDUSTFilteredKmerizer.stream(sb,kSize,maxDUSTScore,kmerSeed)
                            .map(kmer -> kmer.canonical(kSize)).forEach(kmers::add);
                    sb.setLength(0);
                }
            }
            if ( sb.length() > 0 ) {
                SVDUSTFilteredKmerizer.stream(sb,kSize,maxDUSTScore,kmerSeed)
                        .map(kmer -> kmer.canonical(kSize)).forEach(kmers::add);
            }
            return kmers;
        }
        catch ( IOException ioe ) {
            throw new GATKException("Can't read high copy kmers fasta file "+fastaFilename, ioe);
        }
    }

}
