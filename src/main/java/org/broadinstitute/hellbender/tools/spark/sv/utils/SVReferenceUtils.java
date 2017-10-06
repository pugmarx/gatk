package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class SVReferenceUtils {

    public static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    public static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

    public static List<String> getRefNames(final SAMFileHeader header) {
        return header.getSequenceDictionary().getSequences().stream()
                .map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList());
    }

    /**
     * Create an RDD from the reference sequences.
     * The reference sequences are transformed into a single, large collection of byte arrays. The collection is then
     * parallelized into an RDD.
     * Each contig that exceeds a size given by REF_RECORD_LEN is broken into a series of REF_RECORD_LEN chunks with a
     * K-1 base overlap between successive chunks. (I.e., for K=63, the last 62 bases in chunk n match the first 62
     * bases in chunk n+1) so that we don't miss any kmers due to the chunking -- we can just kmerize each record
     * independently.
     */
    public static JavaRDD<byte[]> getRefRDD(final JavaSparkContext ctx,
                                            final int kSize,
                                            final ReferenceMultiSource ref,
                                            final SAMSequenceDictionary readsDict,
                                            final int ref_record_len,
                                            final int ref_records_per_partition) {
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
        if ( dict == null ) throw new GATKException("No reference dictionary available");

        final int effectiveRecLen = ref_record_len - kSize + 1;
        final List<byte[]> sequenceChunks = new ArrayList<>();
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String seqName = rec.getSequenceName();
            final int seqLen = rec.getSequenceLength();
            final SimpleInterval interval = new SimpleInterval(seqName, 1, seqLen);
            try {
                final byte[] bases = ref.getReferenceBases(null, interval).getBases();
                for ( int start = 0; start < seqLen; start += effectiveRecLen ) {
                    sequenceChunks.add(Arrays.copyOfRange(bases, start, Math.min(start+ref_record_len, seqLen)));
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't get reference sequence bases for " + interval, ioe);
            }
        }

        return ctx.parallelize(sequenceChunks, sequenceChunks.size()/ref_records_per_partition+1);
    }

    /**
     * Do a map/reduce on an RDD of genomic sequences:
     * Kmerize, mapping to a pair <kmer,1>, reduce by summing values by key, filter out <kmer,N> where
     * N <= MAX_KMER_FREQ, and collect the high frequency kmers back in the driver.
     */
    @VisibleForTesting
    public static List<SVKmer> processRefRDD(final int kSize,
                                             final int maxDUSTScore,
                                             final int maxKmerFreq,
                                             final JavaRDD<byte[]> refRDD) {
        final int nPartitions = refRDD.getNumPartitions();
        final int hashSize = 2*REF_RECORDS_PER_PARTITION;
        return refRDD
                .mapPartitions(seqItr -> {
                    final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerCounts = new HopscotchMap<>(hashSize);
                    while ( seqItr.hasNext() ) {
                        final byte[] seq = seqItr.next();
                        SVDUSTFilteredKmerizer.canonicalStream(seq, kSize, maxDUSTScore, new SVKmerLong())
                                .forEach(kmer -> {
                                    final KmerAndCount entry = kmerCounts.find(kmer);
                                    if ( entry == null ) kmerCounts.add(new KmerAndCount((SVKmerLong)kmer));
                                    else entry.bumpCount();
                                });
                    }
                    return kmerCounts.iterator();
                })
                .mapToPair(entry -> new Tuple2<>(entry.getKey(), entry.getValue()))
                .partitionBy(new HashPartitioner(nPartitions))
                .mapPartitions(pairItr -> {
                    final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerCounts =
                            new HopscotchMap<>(hashSize);
                    while ( pairItr.hasNext() ) {
                        final Tuple2<SVKmer, Integer> pair = pairItr.next();
                        final SVKmer kmer = pair._1();
                        final int count = pair._2();
                        KmerAndCount entry = kmerCounts.find(kmer);
                        if ( entry == null ) kmerCounts.add(new KmerAndCount((SVKmerLong)kmer, count));
                        else entry.bumpCount(count);
                    }
                    return kmerCounts.stream()
                            .filter(kmerAndCount -> kmerAndCount.grabCount() > maxKmerFreq)
                            .map(KmerAndCount::getKey).iterator();
                })
                .collect();
    }
}
