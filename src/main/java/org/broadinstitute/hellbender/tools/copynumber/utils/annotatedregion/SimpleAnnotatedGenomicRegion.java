package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Map;

/**
 * Simple class that just has an interval and name value pairs.
 */
public class SimpleAnnotatedGenomicRegion implements Locatable {
    private SimpleInterval interval;
    private Map<String, String> annotations;

    public SimpleAnnotatedGenomicRegion(SimpleInterval interval, Map<String, String> annotations) {
        this.interval = interval;
        this.annotations = annotations;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public Map<String, String> getAnnotations() {
        return annotations;
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
}
