package bdsky;

import bdsky.evolution.speciation.JointTreeAndRho;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;

import java.io.PrintStream;

public class MarginalLikelihoodVarianceLogger extends BEASTObject implements Loggable {

    public Input<JointTreeAndRho> modelInput = new Input<>(
            "model",
            "The tree prior for which we estimate the marginal likelihood.",
            Input.Validate.REQUIRED);

    JointTreeAndRho model;

    @Override
    public void initAndValidate() {
        model = modelInput.get();
    }

    @Override
    public void init(PrintStream out) {
        out.print("marginalLikelihoodVariance\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(String.format("%f", model.getMLVariance()) + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
