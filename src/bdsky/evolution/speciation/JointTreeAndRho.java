package bdsky.evolution.speciation;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Extends the sampled-ancestor birth-death skyline model to condition on " +
        "the times of sequentially sampled tips.")

public class JointTreeAndRho extends BirthDeathSkylineModel {
    @Override
    public void initAndValidate() {

        StringBuilder leafTimes = new StringBuilder();

        for (Node leaf : treeInput.get().getExternalNodes()) {
            leafTimes.append(leaf.getHeight()).append(" ");
        }

        rhoSamplingTimes.setValue(leafTimes.toString(), this);// sets the leaf times to the rho sampling times
        reverseTimeArraysInput.setValue("true true true true true",this);// sets reverseTimeArraysInput to true for all time arrays
        conditionOnSurvival.setValue("false",this);

        super.initAndValidate();
    }

    public double calculateSurvivalProbability() {
//        computes the probability that the total population dies out before the final sample time
        double[] ai = new double[totalIntervals];
        double[] bi = new double[totalIntervals];
        double[] P0 = new double[totalIntervals];

        for(int i=0 ; i<totalIntervals ; i++)
            ai[i] = Ai(birth[i],death[i],psi[i]);


        bi[totalIntervals - 1] = Bi(
                birth[totalIntervals - 1],
                death[totalIntervals - 1],
                psi[totalIntervals - 1],
                1.0,
                ai[totalIntervals - 1], 1.0);

        for (int i = totalIntervals - 2; i >= 0; i--) {
            P0[i + 1] = p0(birth[i + 1], death[i + 1], psi[i + 1], ai[i + 1], bi[i+1], times[i + 1], times[i]);

            if (Math.abs(P0[i + 1] - 1) < 1e-10) {
                return Double.NEGATIVE_INFINITY;
            }

            bi[i] = Bi(birth[i], death[i], psi[i], 0.0, ai[i], P0[i + 1]);
        }

        P0[0] = p0(birth[1], death[1], psi[1], ai[1], bi[1], times[1], 0);

        double survivalProb;
        survivalProb = 1 - P0[0];
        return survivalProb;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        logP = super.calculateTreeLogLikelihood(tree); //computes density of the tree given rho

        for (Double thisRho : m_rho.get().getValues()) { // adjusts for the joint density of the tree and rho
            logP = logP - Math.log(thisRho);
        }

        if(SAModel)
            logP = logP - Math.log(calculateSurvivalProbability()); // conditions on the total population size being at least one

        return logP;
    }

    public static void main(String[] args) {

        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:1,B:0.5):1,C:0.5):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        System.out.println(tree);

        JointTreeAndRho jtar = new JointTreeAndRho();
        jtar.initByName(
                "tree", tree,
                "origin", new RealParameter("5.0"),
                "reproductiveNumber", new RealParameter("1.5 1.5"),
                "becomeUninfectiousRate", new RealParameter("1.0"),
                "samplingProportion", new RealParameter("0.0"),
                "rho", new RealParameter("0.1 0.1 0.1"),
                "removalProbability", new RealParameter("0.0"),
                "birthRateChangeTimes", new RealParameter("0.0 2.5")
        );
        System.out.println(jtar.calculateLogP());
        System.out.println(jtar.calculateSurvivalProbability());
    }
}