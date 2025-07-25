package bdsky.util;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;
import beast.base.util.Randomizer;

import java.util.HashSet;
import java.util.Set;

@Description("Extension of TraitSet which ensures trait values are unique" +
        " by adding a small random perturbation to each.")
public class PerturbedTraitSet extends TraitSet {

    public Input<Double> epsilonInput = new Input<>("epsilon",
            "Maximum perturbation to apply to trait values.",
            1e-5);

    public void initAndValidate() {
        super.initAndValidate();

        double epsilon = epsilonInput.get();

        Set<Double> seenValues = new HashSet<>();
        for (int i=0; i<values.length; i++) {
            if (seenValues.contains(values[i])) {
                double perturbedValue;
                do {
                    perturbedValue = values[i]
                            - 0.5*epsilon + Randomizer.nextDouble()*epsilon;
                } while (seenValues.contains(perturbedValue));
                values[i] = perturbedValue;
            } else {
                seenValues.add(values[i]);
            }
        }

        // update extremes
        minValue = values[0];
        maxValue = values[0];
        for (double value : values) {
            minValue = Math.min(minValue, value);
            maxValue = Math.max(maxValue, value);
        }
    }
}
