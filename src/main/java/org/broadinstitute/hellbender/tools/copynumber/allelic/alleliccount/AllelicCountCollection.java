package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LocatableCollectionTSV;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Simple data structure to pass and read/write a List of {@link AllelicCount} objects.
 * All {@link AllelicCount} fields (including ref/alt nucleotide) must be specified if reading/writing from/to file.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class AllelicCountCollection extends LocatableCollectionTSV<AllelicCount> {
    enum AllelicCountTableColumn {
        CONTIG,
        POSITION,
        REF_COUNT,
        ALT_COUNT,
        REF_NUCLEOTIDE,
        ALT_NUCLEOTIDE;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
    
    private static final Function<DataLine, AllelicCount> ALLELIC_COUNT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(AllelicCountTableColumn.CONTIG);
        final int position = dataLine.getInt(AllelicCountTableColumn.POSITION);
        final int refReadCount = dataLine.getInt(AllelicCountTableColumn.REF_COUNT);
        final int altReadCount = dataLine.getInt(AllelicCountTableColumn.ALT_COUNT);
        final Nucleotide refNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumn.REF_NUCLEOTIDE.name()).getBytes()[0]);
        final Nucleotide altNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumn.ALT_NUCLEOTIDE.name()).getBytes()[0]);
        final SimpleInterval interval = new SimpleInterval(contig, position, position);
        return new AllelicCount(interval, refReadCount, altReadCount, refNucleotide, altNucleotide);
    };

    private static final BiConsumer<AllelicCount, DataLine> ALLELIC_COUNT_RECORD_AND_DATA_LINE_BI_CONSUMER = (allelicCount, dataLine) ->
            dataLine.append(allelicCount.getInterval().getContig())
                    .append(allelicCount.getInterval().getEnd())
                    .append(allelicCount.getRefReadCount())
                    .append(allelicCount.getAltReadCount())
                    .append(allelicCount.getRefNucleotide().name())
                    .append(allelicCount.getAltNucleotide().name());

    public AllelicCountCollection(final File inputFile) {
        super(inputFile, AllelicCountCollection.AllelicCountTableColumn.COLUMNS, ALLELIC_COUNT_DATA_LINE_TO_RECORD_FUNCTION, ALLELIC_COUNT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public AllelicCountCollection(final String sampleName,
                                  final List<AllelicCount> AllelicCounts) {
        super(sampleName, AllelicCounts, AllelicCountCollection.AllelicCountTableColumn.COLUMNS, ALLELIC_COUNT_DATA_LINE_TO_RECORD_FUNCTION, ALLELIC_COUNT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }
}
