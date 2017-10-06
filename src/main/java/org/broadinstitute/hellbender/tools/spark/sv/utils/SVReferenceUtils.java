package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import scala.Tuple2;

import java.util.List;
import java.util.stream.Collectors;

public final class SVReferenceUtils {
    public static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    public static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

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

    public static List<String> getRefNames(final SAMFileHeader header) {
        return header.getSequenceDictionary().getSequences().stream()
                .map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList());
    }
}
