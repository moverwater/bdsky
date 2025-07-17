package bdsky.evolution.speciation;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.RealParameter;

import java.util.Comparator;
import java.util.List;

public class TreeTipPerturber extends BEASTObject implements StateNodeInitialiser {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree to initialize",
            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {

    }

    @Override
    public void initStateNodes() {
        Tree tree = treeInput.get();

//      Perturb tip heights
        List<Node> leaves = tree.getExternalNodes();
        double epsilon = 1e-4;

// Sort leaves by height (smallest first)
        leaves.sort(Comparator.comparingDouble(Node::getHeight));

        for (int i = 0; i < leaves.size(); i++) {
            Node leaf = leaves.get(i);
            double originalHeight = leaf.getHeight();
//            double randomOffset = (2 * Math.random() - 1) * epsilon;
            double newHeight = originalHeight + i * epsilon;
            leaf.setHeight(newHeight); // Apply small offset
        }
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {

    }
}
