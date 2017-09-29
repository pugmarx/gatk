package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A class for storing end-end integration tests over artificial data intended to test allele specific annotation implementations,
 * As of 8/29/17, test output was matched against GATK3 implementations of combineGVCFs().
 */
public abstract class ReducibleAnnotationBaseTest extends BaseTest {

    @Override
    public String getToolTestDataDir() {
        return "src/test/resources/" + this.getClass().getPackage().getName().replace(".", "/") + "/";
    }

    @DataProvider(name = "interestingSitesCombineResults")
    public Object[][] makeAnnotateRsIDData() {
        FeatureDataSource<VariantContext> VCFA = new FeatureDataSource<>(getTestFile("NA12878.AS.chr20snippet.g.vcf"));
        FeatureDataSource<VariantContext> VCFB = new FeatureDataSource<>(getTestFile("NA12892.AS.chr20snippet.g.vcf"));
        FeatureDataSource<VariantContext> CombineVCFOutput = new FeatureDataSource<>(getTestFile("CombineGVCFs.output.vcf"));
        FeatureDataSource<VariantContext> GenotypeGVCFOutput = new FeatureDataSource<>(getTestFile("GenotypeGVCFs.output.vcf"));

        List<Object[]> tests = new ArrayList<>();

        // these are hand picked sites from the allele specific unit tests for combinegvcfs that triggered a combine in GATK3
        Integer[] interestingLocs = {10087820, 10433312, 10433322, 10433324, 10433326, 10433328, 10433382, 10433391, 10433468, 10433560, 10433575, 10433594, 10433955,
                10434177, 10434384, 10435067, 10434258, 10436227, 10684106};
        List<SimpleInterval> intervals = Arrays.stream(interestingLocs).map(m -> new SimpleInterval("20", m, m)).collect(Collectors.toList());
        for (SimpleInterval loc : intervals) {
            VariantContext a = VCFA.query(loc).next();
            VariantContext b = VCFB.query(loc).next();
            VariantContext result = CombineVCFOutput.query(loc).next();
            Iterator<VariantContext> query = GenotypeGVCFOutput.query(loc);
            VariantContext genotyped = query.hasNext() ? query.next() : null;
            tests.add(new Object[]{Arrays.asList(a, b), result, genotyped});
        }
        VCFA.close();
        VCFB.close();
        CombineVCFOutput.close();
        GenotypeGVCFOutput.close();

        return tests.toArray(new Object[][]{});
    }

    /*
     * Methods that must be overridden in order for the automatic GATK3 combineGVCFs tests to be run.
     */
    protected abstract List<String> getAnnotationsToUse();

    protected abstract String getRawKey();

    protected abstract String getKey();


    // NOTE: this code is mimicking the behavior of GATK3 combineGVCFS insofar as it is important for the annotations
    @Test(dataProvider = "interestingSitesCombineResults")
    public void testCombineAnnotationGATK3Concordance(List<VariantContext> VCs, VariantContext result, VariantContext genotyped) throws Exception {
        VariantAnnotatorEngine annotatorEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(Collections.emptyList(), getAnnotationsToUse(), Collections.emptyList(), null, Collections.emptyList());
        ReferenceConfidenceVariantContextMerger merger = new ReferenceConfidenceVariantContextMerger(annotatorEngine);
        VariantContext merged = merger.merge(VCs, new SimpleInterval(result.getContig(), result.getStart(), result.getStart()), result.getReference().getBases()[0], false, false);
        Assert.assertTrue(VariantContextTestUtils.alleleSpecificAnnotationEquals(merged, result, getRawKey()));
    }

    // NOTE: this code is mimicking the behavior of GATK3 GenotypeGVCFs
    @Test(dataProvider = "interestingSitesCombineResults")
    public void testFinalizeAnnotationGATK3Concordance(List<VariantContext> VCs, VariantContext result, VariantContext genotyped) throws Exception {
        if (result == null || genotyped == null) {
            return;
        }

        VariantAnnotatorEngine annotatorEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(Collections.emptyList(), getAnnotationsToUse(), Collections.emptyList(), null, Collections.emptyList());
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = new GenotypeCalculationArgumentCollection();
        GeneralPloidyFailOverAFCalculatorProvider calculatorProvider = new GeneralPloidyFailOverAFCalculatorProvider(uac.genotypeArgs);

        GenotypingEngine<?> genotypingEngine = new MinimalGenotypingEngine(uac, new IndexedSampleList(result.getSampleNamesOrderedByName()), calculatorProvider);
        genotypingEngine.setAnnotationEngine(annotatorEngine);
        GenotypeLikelihoodsCalculationModel model = result.getType() == VariantContext.Type.INDEL
                ? GenotypeLikelihoodsCalculationModel.INDEL
                : GenotypeLikelihoodsCalculationModel.SNP;
        VariantContext withGenotypes = genotypingEngine.calculateGenotypes(result, model, null);
        withGenotypes = new VariantContextBuilder(withGenotypes).attributes(result.getAttributes()).make();
        VariantContext finalized = annotatorEngine.finalizeAnnotations(withGenotypes, result);
        finalized =  GATKVariantContextUtils.reverseTrimAlleles(finalized);
        Assert.assertTrue(VariantContextTestUtils.alleleSpecificAnnotationEquals(finalized, genotyped, getKey()));
    }
}
