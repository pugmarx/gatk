package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;

import java.math.BigInteger;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * Useful scraps of this and that.
 */
public final class SVUtils {

    // =================================================================================================================

    /** return a good initialCapacity for a HashMap that will hold a given number of elements */
    public static int hashMapCapacity( final int nElements )
    {
        return (int)((nElements*4L)/3) + 1;
    }

    /** count the number of items available from an iterator */
    public static <T> int iteratorSize( final Iterator<T> itr ) {
        int result = 0;
        while ( itr.hasNext() ) { result += 1; itr.next(); }
        return result;
    }

    public static int iteratorSize( final LongIterator itr ) {
        int result = 0;
        while ( itr.hasNext() ) { result += 1; itr.next(); }
        return result;
    }

    public static <T> Iterator<T> singletonIterator( final T t ) {
        return Collections.singletonList(t).iterator();
    }

    public static Collection<SVKmer> uniquify(final Collection<SVKmer> coll1, final Collection<SVKmer> coll2) {
        final HopscotchSet<SVKmer> kmers = new HopscotchSet<>(coll1.size() + coll2.size());
        kmers.addAll(coll1);
        kmers.addAll(coll2);
        return kmers;
    }

    public static class IteratorFilter<T> implements Iterator<T> {
        private final Iterator<T> itr;
        private final Predicate<T> predicate;
        private T obj;

        public IteratorFilter( final Iterator<T> itr, final Predicate<T> predicate ) {
            this.itr = itr;
            this.predicate = predicate;
            advance();
        }

        @Override public boolean hasNext() { return obj != null; }

        @Override
        public T next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException("IteratorFilter is exhausted.");
            }
            final T result = obj;
            advance();
            return result;
        }

        private void advance() {
            obj = null;
            while ( itr.hasNext() ) {
                final T next = itr.next();
                if ( predicate.test(next) ) {
                    obj = next;
                    break;
                }
            }
        }
    }

    /**
     * Provides a stream collector that will collect items into an array list with a given initial capacity.
     */
    public static <T> Collector<T, ?, ArrayList<T>> arrayListCollector(final int size) {
        return Collectors.toCollection( () -> new ArrayList<>(size));
    }

    // =================================================================================================================

    //Workaround for seed 14695981039346656037 that doesn't fit in a signed long
    private static final long FNV64_DEFAULT_SEED = new BigInteger("14695981039346656037").longValue();

    /**
     * 64-bit FNV-1a hash for long's
     */

    public static long fnvLong64( final long toHash ) {
        return fnvLong64(FNV64_DEFAULT_SEED, toHash);
    }

    public static long fnvLong64( long start, final long toHash ) {
        final long mult = 1099511628211L;
        start ^= (toHash >> 56) & 0xffL;
        start *= mult;
        start ^= (toHash >> 48) & 0xffL;
        start *= mult;
        start ^= (toHash >> 40) & 0xffL;
        start *= mult;
        start ^= (toHash >> 32) & 0xffL;
        start *= mult;
        start ^= (toHash >> 24) & 0xffL;
        start *= mult;
        start ^= (toHash >> 16) & 0xffL;
        start *= mult;
        start ^= (toHash >> 8) & 0xffL;
        start *= mult;
        start ^= toHash & 0xffL;
        start *= mult;
        return start;
    }

    /**
     * 64-bit FNV-1a hash for byte arrays
     */
    public static long fnvByteArray64(final byte[] toHash) {
        // TODO: this is a mistake:  the constant should be the FNV64_DEFAULT_SEED, but it's the multiplier instead.
        return fnvByteArray64(1099511628211L, toHash);
    }

    public static long fnvByteArray64(long start, final byte[] toHash) {
        for (int i = 0; i < toHash.length; i += 8) {
            long val = 0;
            for (int j = 0; j < 8 && i + j < toHash.length; j++) {
                val = (val << 8) | toHash[i + j];
            }
            start = fnvLong64(start, val);
        }
        return start;
    }

}
