package org.broadinstitute.hellbender.tools.copynumber.legacy.formats;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Unit tests for {@link LocatableCollectionTSV}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class LocatableCollectionTSVUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/legacy/formats";
    private static final File SIMPLE_LOCATABLE_COLLECTION_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection.tsv");
    private static final File SIMPLE_LOCATABLE_COLLECTION_NON_LEXICOGRAPHICAL_ORDER_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection-non-lexicographical-order.tsv");
    private static final String SAMPLE_NAME_EXPECTED = "test";

    //simple example of a record class
    private static final class SimpleLocatable implements Locatable {
        private final SimpleInterval interval;
        private final int value;

        private SimpleLocatable(final SimpleInterval interval, final int value) {
            this.interval = Utils.nonNull(interval);
            this.value = value;
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }

        public SimpleInterval getInterval() {
            return interval;
        }

        public double getValue() {
            return value;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final SimpleLocatable that = (SimpleLocatable) o;

            return value == that.value && interval.equals(that.interval);
        }

        @Override
        public int hashCode() {
            int result = interval.hashCode();
            result = 31 * result + value;
            return result;
        }
    }

    //simple example of a collection class
    private static final class SimpleLocatableCollection extends LocatableCollectionTSV<SimpleLocatable> {
        enum SimpleLocatableTableColumn {
            CONTIG,
            START,
            END,
            VALUE;

            static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
        }
        
        private static final Function<DataLine, SimpleLocatable> SIMPLE_LOCATABLE_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
            final String contig = dataLine.get(SimpleLocatableTableColumn.CONTIG);
            final int start = dataLine.getInt(SimpleLocatableTableColumn.START);
            final int end = dataLine.getInt(SimpleLocatableTableColumn.END);
            final int value = dataLine.getInt(SimpleLocatableTableColumn.VALUE);
            final SimpleInterval interval = new SimpleInterval(contig, start, end);
            return new SimpleLocatable(interval, value);
        };

        private static final BiConsumer<SimpleLocatable, DataLine> SIMPLE_LOCATABLE_RECORD_AND_DATA_LINE_BI_CONSUMER = (simpleLocatable, dataLine) ->
                dataLine.append(simpleLocatable.getInterval().getContig())
                        .append(simpleLocatable.getInterval().getStart())
                        .append(simpleLocatable.getInterval().getEnd())
                        .append(simpleLocatable.getValue());

        private SimpleLocatableCollection(final File inputFile) {
            super(inputFile, SimpleLocatableTableColumn.COLUMNS, SIMPLE_LOCATABLE_DATA_LINE_TO_RECORD_FUNCTION, SIMPLE_LOCATABLE_RECORD_AND_DATA_LINE_BI_CONSUMER);
        }

        private SimpleLocatableCollection(final String sampleName,
                                         final List<SimpleLocatable> simpleLocatables) {
            super(sampleName, simpleLocatables, SimpleLocatableTableColumn.COLUMNS, SIMPLE_LOCATABLE_DATA_LINE_TO_RECORD_FUNCTION, SIMPLE_LOCATABLE_RECORD_AND_DATA_LINE_BI_CONSUMER);
        }
    }

    private static final SimpleLocatableCollection SIMPLE_LOCATABLE_COLLECTION_EXPECTED = new SimpleLocatableCollection(
            SAMPLE_NAME_EXPECTED,
            Arrays.asList(
                    new SimpleLocatable(new SimpleInterval("1", 1, 1), 1),
                    new SimpleLocatable(new SimpleInterval("1", 2, 2), 2),
                    new SimpleLocatable(new SimpleInterval("10", 1, 1), 3),
                    new SimpleLocatable(new SimpleInterval("2", 1, 1), 4)));

    @Test
    public void testRead() {
        final SimpleLocatableCollection simpleLocatableCollection = new SimpleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_FILE);
        Assert.assertFalse(simpleLocatableCollection == SIMPLE_LOCATABLE_COLLECTION_EXPECTED);
        Assert.assertEquals(simpleLocatableCollection, SIMPLE_LOCATABLE_COLLECTION_EXPECTED);
    }

    @Test
    public void testWrite() throws IOException {
        final File tempFile = createTempFile("test", ".tsv");
        SIMPLE_LOCATABLE_COLLECTION_EXPECTED.write(tempFile);
        Assert.assertTrue(FileUtils.contentEquals(tempFile, SIMPLE_LOCATABLE_COLLECTION_FILE));
    }

    @Test
    public void testLexicographicalSortingByConstructor() {
        final SimpleLocatableCollection simpleLocatableCollectionExpectedUnsortedListArgument = new SimpleLocatableCollection(
                SAMPLE_NAME_EXPECTED,
                Arrays.asList(
                        new SimpleLocatable(new SimpleInterval("1", 1, 1), 1),
                        new SimpleLocatable(new SimpleInterval("1", 2, 2), 2),
                        new SimpleLocatable(new SimpleInterval("2", 1, 1), 4),
                        new SimpleLocatable(new SimpleInterval("10", 1, 1), 3)));
        Assert.assertFalse(simpleLocatableCollectionExpectedUnsortedListArgument == SIMPLE_LOCATABLE_COLLECTION_EXPECTED);
        Assert.assertEquals(simpleLocatableCollectionExpectedUnsortedListArgument, SIMPLE_LOCATABLE_COLLECTION_EXPECTED);
    }
    
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBadOrderInput() {
        new SimpleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_NON_LEXICOGRAPHICAL_ORDER_FILE);
    }
}