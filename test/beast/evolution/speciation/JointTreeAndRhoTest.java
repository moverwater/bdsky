package beast.evolution.speciation;

import bdsky.evolution.speciation.JointTreeAndRho;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

public class JointTreeAndRhoTest {

    @Test
    public void test() {
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:1,B:2):1,C:1):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        System.out.println(tree);

        JointTreeAndRho jtar = new JointTreeAndRho();
        jtar.initByName(
                "tree", tree,
                "origin", new RealParameter("5.0"),
                "reproductiveNumber", new RealParameter("0.2"),
                "becomeUninfectiousRate", new RealParameter("1.0"),
                "samplingProportion", new RealParameter("0.0"),
                "rho", new RealParameter("0.1 0.1 0.1"),
                "removalProbability", new RealParameter("0.0")
        );
//        System.out.println(jtar.totalIntervals);
        System.out.println("logP=" + jtar.calculateLogP());
        System.out.println("1-p0=" + jtar.calculateSurvivalProbability());
        System.out.println("leafTimes" + Arrays.toString(jtar.rhoSamplingTimes.get().getValues()));

        assertEquals(1.3, jtar.calculateLogP(), 1e-10);
    }
}
