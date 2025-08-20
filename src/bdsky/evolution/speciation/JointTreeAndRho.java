package bdsky.evolution.speciation;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * @author Marcus Overwater
 */

@Description("Extends the sampled-ancestor birth-death skyline model to condition on " +
        "the times of sequentially sampled tips.")

public class JointTreeAndRho extends BirthDeathSkylineModel {
    public Input<Boolean> approximateMarginalLikelihoodInput =
            new Input<Boolean>("approximateMarginalLikelihoodInput", "Boolean, true if you want the approximate the marginal likelihood with a monte-carlo estimate, otherwise sample rho as part of the state space", false);

    public Input<IntegerParameter> numberOfMonteCarloSamplesInput =
            new Input<IntegerParameter>("numberOfMonteCarloSamples","The number of samples in the Monte Carlo estimate of the marginal likelihood if approximateMarginalLikelihoodInput is true, default 1000", (IntegerParameter) null);

    protected Double mlVar;
    private Boolean approxMarginal;
    private Integer numSamples;
    protected TreeInterface tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

// Sort leaves by height (smallest first)
        double[] leafHeights = tree.getExternalNodes().stream()
                .mapToDouble(Node::getHeight).sorted().toArray();

        StringBuilder leafTimeString = new StringBuilder();
        for (double leafHeight : leafHeights)
            leafTimeString.append(leafHeight).append(" ");
        rhoSamplingTimes.setValue(leafTimeString.toString(),this);

        approxMarginal = approximateMarginalLikelihoodInput.get();

        if (approxMarginal) {
            if (numberOfMonteCarloSamplesInput.get() != null) {
                numSamples = numberOfMonteCarloSamplesInput.get().getValue();
            }
            else {
                numSamples = 1000;
            }
        }

        reverseTimeArraysInput.setValue("true true true true true",this);// sets reverseTimeArraysInput to true for all time arrays
        conditionOnSurvival.setValue("false",this);
        super.initAndValidate();
    }

    public double calculateSurvivalProbability() {
//        computes the probability that the total population dies out before the final sample time
        double[] ai = new double[totalIntervals];
        double[] bi = new double[totalIntervals];
        double[] P0 = new double[totalIntervals];
            for (int i = 0; i < totalIntervals; i++)
                ai[i] = Ai(birth[i], death[i], psi[i]);

            bi[totalIntervals - 1] = Bi(
                    birth[totalIntervals - 1],
                    death[totalIntervals - 1],
                    psi[totalIntervals - 1],
                    1.0,
                    ai[totalIntervals - 1], 1.0);

            for (int i = totalIntervals - 2; i >= 0; i--) {
                P0[i + 1] = p0(birth[i + 1], death[i + 1], psi[i + 1], ai[i + 1], bi[i + 1], times[i + 1], times[i]);

//                if (Math.abs(P0[i + 1] - 1) < 1e-10) {
//                    return Double.NEGATIVE_INFINITY;
//                }

                bi[i] = Bi(birth[i], death[i], psi[i], 0.0, ai[i], P0[i + 1]);
            }

            P0[0] = p0(birth[1], death[1], psi[1], ai[1], bi[1], times[1], 0);
        double survivalProb;
        survivalProb = 1 - P0[0];
        if (survivalProb == 0.0){
            return Double.NEGATIVE_INFINITY;
        }
        return survivalProb;
    }

    public double calculateJointTreeAndRho(TreeInterface tree) {
        logP = super.calculateTreeLogLikelihood(tree); //computes density of the tree given rho

        for (Double thisRho : m_rho.get().getDoubleValues()) { // adjusts for the joint density of the tree and rho
            logP = logP - Math.log(thisRho);
        }
        if (SAModel)
            logP = logP - Math.log(calculateSurvivalProbability()); // conditions on the total population size being at least one
        return logP;
    }

    @Override
    public boolean isStochastic() {
        return approxMarginal;
    }

    public double getMLVariance() {
        return mlVar;
    }

//    int count=0;

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        if (approxMarginal) {

            int dim = rhoSamplingTimes.get().getDimension();

            double[] rhoSample = new double[dim];

            double maxLogP = Double.NEGATIVE_INFINITY;
            double[] logPs = new double[numSamples];


            int validSamples = 0;

            for (int i = 0; i < numSamples; i++) {
                for (int d = 0; d < dim; d++) {
                    rhoSample[d] = Math.random(); // uniform [0,1]
                }

                StringBuilder sb = new StringBuilder();
                for (int d = 0; d < dim; d++) {
                    sb.append(rhoSample[d]);
                    if (d < dim - 1) {
                        sb.append(" ");
                    }
                }
                m_rho.setValue(sb.toString(), this);

                double logP = calculateJointTreeAndRho(tree);

                if (!Double.isInfinite(logP) && !Double.isNaN(logP)) {
                    logPs[validSamples++] = logP;
                    if (logP > maxLogP) {
                        maxLogP = logP;
                    }
                }
            }

            if (validSamples == 0) {
                return Double.NEGATIVE_INFINITY;
            }

            double sumExp = 0.0;
            double sumLogW = 0.0;
            double sumLogWSq = 0.0;

            for (int i = 0; i < validSamples; i++) {
                double logW = logPs[i] - maxLogP;  // log w_i
                double w = Math.exp(logW);         // exp(log w_i)

                sumExp += w;                       // For logZ
                sumLogW += logW;                   // For mean(log w_i)
                sumLogWSq += logW * logW;          // For variance(log w_i)
            }

            logP = Math.log(sumExp / validSamples) + maxLogP;

// Variance of log w_i
            double meanLogW = sumLogW / validSamples;
            double varLogW = (sumLogWSq - validSamples * meanLogW * meanLogW) / (validSamples - 1);

// Variance of log Z estimate via Delta method in log-space
            mlVar = varLogW / validSamples;
        }
        else {
            logP = calculateJointTreeAndRho(tree);
        }
        return logP;
    }

    public static void main(String[] args) {

        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:2,B:2):1,C:3):0;",
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
        System.out.println(jtar.totalIntervals);
        System.out.println("logP=" + jtar.calculateLogP());
        System.out.println("1-p0=" + jtar.calculateSurvivalProbability());
        System.out.println("isRhoTip=" + Arrays.toString(jtar.isRhoTip));
    }
}