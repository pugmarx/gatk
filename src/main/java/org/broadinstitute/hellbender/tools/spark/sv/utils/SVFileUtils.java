package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.*;
import java.util.*;

public final class SVFileUtils {

    private static final String REFERENCE_GAP_INTERVAL_FILE_COMMENT_LINE_PROMPT = "#";

    /**
     * Creates a directory, in local FS, HDFS, or Google buckets to write to.
     */
    public static boolean createDirToWriteTo(final String pathString) {
        try {
            Utils.nonNull(pathString);
            if ( java.nio.file.Files.exists(java.nio.file.Paths.get(pathString)) )
                throw new IOException("Directory to be created already exists: " + pathString);

            final boolean isSuccessful;
            if (BucketUtils.isHadoopUrl(pathString)) {
                isSuccessful = org.apache.hadoop.fs.FileSystem.get(new org.apache.hadoop.conf.Configuration()).mkdirs(new org.apache.hadoop.fs.Path(pathString));
            } else {
                final java.nio.file.Path dir = java.nio.file.Files.createDirectory(IOUtils.getPath(pathString));
                isSuccessful = java.nio.file.Files.isDirectory(dir) && java.nio.file.Files.isWritable(dir);
            }
            return isSuccessful;
        } catch (final IOException x) {
            throw new GATKException("Could not create directory: " + pathString, x);
        }
    }

    public static void writeLinesToSingleFile(final Iterator<String> linesToWrite, final String fileName) {
        try ( final OutputStream writer =
                      new BufferedOutputStream(BucketUtils.createFile(fileName)) ) {
            while (linesToWrite.hasNext()) {
                writer.write(linesToWrite.next().getBytes()); writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }

    public static void writeSAMFile(final Iterator<SAMRecord> alignments, final SAMFileHeader header, final String outputName,
                                    final boolean preOrdered) {
        try (final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(new File(outputName), null, header,
                preOrdered, true, false) ) {
            alignments.forEachRemaining(writer::addAlignment);
        } catch ( final UncheckedIOException ie) {
            throw new GATKException("Can't write SAM file to the specified location: " + outputName, ie);
        }
    }


    /**
     * Read a file of kmers.
     * Each line must be exactly
     * {@link org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection#KMER_SIZE}
     * characters long, and must match [ACGT]*.
     */
    public static Set<SVKmer> readKmersFile(final int kSize, final String kmersFile,
                                            final SVKmer kmer ) {
        final Set<SVKmer> kmers;

        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(kmersFile))) ) {
            final long fileLength = BucketUtils.fileSize(kmersFile);
            kmers = new HopscotchSet<>((int)(fileLength/(kSize+1)));
            String line;
            while ( (line = rdr.readLine()) != null ) {
                if ( line.length() != kSize ) {
                    throw new GATKException("SVKmer kill set contains a line of length " + line.length() +
                            " but we were expecting K=" + kSize);
                }

                final SVKmerizer kmerizer = new SVKmerizer(line, kSize, 1, new SVKmerLong(kSize));
                if ( !kmerizer.hasNext() ) {
                    throw new GATKException("Unable to kmerize the kmer kill set string '" + line + "'.");
                }

                kmers.add(kmerizer.next());
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read kmers from "+kmersFile, ioe);
        }

        return kmers;
    }

    /** Write kmers to file. */
    public static <KType extends SVKmer> void writeKmersFile(final int kSize, final String kmersFile,
                                                             final Collection<KType> kmers ) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(kmersFile))) ) {
            for ( final KType kmer : kmers ) {
                writer.write(kmer.toString(kSize));
                writer.write('\n');
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to write kmers to "+kmersFile, ioe);
        }
    }

    /** Read intervals from file. */
    public static List<SVInterval> readIntervalsFile(final String intervalsFile,
                                                     final Map<String, Integer> contigNameMap ) {
        final List<SVInterval> intervals;
        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(intervalsFile))) ) {
            final int INTERVAL_FILE_LINE_LENGTH_GUESS = 25;
            final long sizeGuess = BucketUtils.fileSize(intervalsFile)/INTERVAL_FILE_LINE_LENGTH_GUESS;
            intervals = new ArrayList<>((int)sizeGuess);
            String line;
            int lineNo = 0;
            while ( (line = rdr.readLine()) != null ) {
                ++lineNo;
                if (line.startsWith(REFERENCE_GAP_INTERVAL_FILE_COMMENT_LINE_PROMPT)) {
                    continue;
                }
                final String[] tokens = line.split("\t");
                if ( tokens.length != 3 ) {
                    throw new GATKException("Interval file "+intervalsFile+" line "+
                            lineNo+" did not contain 3 columns: "+line);
                }
                try {
                    final Integer contigId = contigNameMap.get(tokens[0]);
                    if ( contigId == null ) throw new GATKException("contig name "+tokens[0]+" not in dictionary");
                    final int start = Integer.valueOf(tokens[1]);
                    final int end = Integer.valueOf(tokens[2]);
                    intervals.add(new SVInterval(contigId, start, end));
                }
                catch ( final Exception e ) {
                    throw new GATKException("Unable to parse interval file "+intervalsFile+" line "+lineNo+": "+line, e);
                }
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read intervals from "+intervalsFile, ioe);
        }
        return intervals;
    }

    /** Write intervals to a file. */
    public static void writeIntervalsFile( final String intervalsFile,
                                           final Collection<SVInterval> intervals, final List<String> contigNames ) {
        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(intervalsFile)))) {
            for (final SVInterval interval : intervals) {
                final String seqName = contigNames.get(interval.getContig());
                writer.write(seqName + "\t" + interval.getStart() + "\t" + interval.getEnd() + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write intervals file " + intervalsFile, ioe);
        }
    }
}
