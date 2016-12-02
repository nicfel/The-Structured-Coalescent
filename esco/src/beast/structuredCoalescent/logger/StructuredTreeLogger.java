package beast.structuredCoalescent.logger;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.StateNode;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;

import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.structuredCoalescent.distribution.ApproximateStructuredCoalescent;
import beast.structuredCoalescent.distribution.ExactStructuredCoalescent;
import beast.structuredCoalescent.distribution.IndependentStructuredCoalescent;


/**
 * adapted by Nicola Felix Mueller from the tree logger
 */
@Description("log trees that also contain the node state probabilities")
public class StructuredTreeLogger extends Tree implements Loggable {
	public Input<ExactStructuredCoalescent> exactInput = new Input<ExactStructuredCoalescent>(
			"exactDensity",
			"A deterministic epidemiological model");
	public Input<IndependentStructuredCoalescent> independentInput = new Input<IndependentStructuredCoalescent>(
			"independentDensity",
			"A deterministic epidemiological model");
	public Input<ApproximateStructuredCoalescent> approximateInput = new Input<ApproximateStructuredCoalescent>(
			"approximateDensity",
			"A deterministic epidemiological model");	
	
    public Input<Tree> treeInput = new Input<Tree>("tree", "tree to be logged", Validate.REQUIRED);

    public Input<BranchRateModel.Base> clockModelInput = new Input<BranchRateModel.Base>("branchratemodel", "rate to be logged with branches of the tree");
    public Input<List<Function>> parameterInput = new Input<List<Function>>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    public Input<Boolean> maxStateInput = new Input<Boolean>("maxState", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<BooleanParameter> conditionalStateProbsInput = new Input<BooleanParameter>("conditionalStateProbs", "report branch lengths as substitutions (branch length times clock rate for the branch)");
    public Input<Boolean> substitutionsInput = new Input<Boolean>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Integer> decimalPlacesInput = new Input<Integer>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);
    public Input<TraitSet> typeTraitInput = new Input<>(
            "typeTrait", "Type trait set.");
    
    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;
    boolean takeMax = true;
    boolean conditionals = true;

    private DecimalFormat df;
    private String type;
    private boolean traitInput = false;
    
    private int states;
	
	
    @Override
    public void initAndValidate() {
    	
    	if (conditionalStateProbsInput.get() != null) conditionals = conditionalStateProbsInput.get().getValue();
    	
        if (typeTraitInput.get() != null) traitInput = true;

    	
        if (parameterInput.get().size() == 0 && clockModelInput.get() == null) {
        	someMetaDataNeedsLogging = false;
        	return;
            //throw new Exception("At least one of the metadata and branchratemodel inputs must be defined");
        }
    	someMetaDataNeedsLogging = true;
    	// without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
        	substitutions = substitutionsInput.get();
        }
       
        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();
        	
        }

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }
        
        if(exactInput.get()!=null)
        	type = exactInput.get().getType();
        else if(independentInput.get()!=null)
        	type = independentInput.get().getType();
        else if(approximateInput.get()!=null)
        	type = approximateInput.get().getType();
        else
        	System.err.println("PROBLEM IN STRUCTURED LOGGER");
        
        
        states = 0;
        
        if(exactInput.get()!=null)
        	states = exactInput.get().states;
        else if(independentInput.get()!=null)
        	states = independentInput.get().states;
        else if(approximateInput.get()!=null)
        	states = approximateInput.get().states;
        
        System.out.println(states);
        System.exit(0);

    }

    @Override
    public void init(PrintStream out) {
    	treeInput.get().init(out);
    }

    @Override
    public void log(int nSample, PrintStream out) {
    	
        states = 0;
        
        if(exactInput.get()!=null)
        	states = exactInput.get().states;
        else if(independentInput.get()!=null)
        	states = independentInput.get().states;
        else if(approximateInput.get()!=null)
        	states = approximateInput.get().states;
        
//        System.out.println(states);
//        System.exit(0);

    	
    	// calculate the true state probabilities to print them
//    	if (conditionals)
//    		volzDensityInput.get().calculateConditionalStates();
//    	else
    	
        if(exactInput.get()!=null)
        	exactInput.get().calculateLogP();
        else if(independentInput.get()!=null)
        	independentInput.get().calculateLogP();
        else if(approximateInput.get()!=null)
        	approximateInput.get().calculateLogP();;

    	
    	
        // make sure we get the current version of the inputs
        Tree tree = (Tree) treeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
//        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadata, branchRateModel));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }

    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();
        	
        }
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
        if (!node.isLeaf()) {
        	if (!takeMax){	        
		        buf.append("[&" + type + "prob={");
		        
		        DoubleMatrix stateProbs = new DoubleMatrix();
		        
		        if(exactInput.get()!=null)
		        	stateProbs = exactInput.get().getStateProb(node.getNr());
		        else if(independentInput.get()!=null)
		        	stateProbs = independentInput.get().getStateProb(node.getNr());
		        else if(approximateInput.get()!=null)
		        	stateProbs = approximateInput.get().getStateProb(node.getNr());

		        
		        
		        for (int i = 0 ; i < states-1; i++){
		        	buf.append(String.format("%.3f", stateProbs.get(i)));
		        	buf.append(",");
		        }
		        buf.append(String.format("%.3f", stateProbs.get(states-1)));
		        buf.append("}");
	        
		        buf.append(",max" + type + "=");
		        buf.append(String.format("%d", stateProbs.argmax() ));
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");
		        DoubleMatrix stateProbs = new DoubleMatrix();
		        
		        if(exactInput.get()!=null)
		        	stateProbs = exactInput.get().getStateProb(node.getNr());
		        else if(independentInput.get()!=null)
		        	stateProbs = independentInput.get().getStateProb(node.getNr());
		        else if(approximateInput.get()!=null)
		        	stateProbs = approximateInput.get().getStateProb(node.getNr());

		        buf.append(String.format("%d", stateProbs.argmax() ));
		        buf.append(']');        		
        	}
        }else{
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;
			
			if(traitInput){				
				sampleState = (int) typeTraitInput.get().getValue(node.getID());
			}
			else{
				sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
			}
			if (!takeMax){
    	        
		        buf.append("[&" + type + "prob={");
	
		        for (int i = 0 ; i < states-1; i++){
		        	if (sampleState != i) buf.append(String.format("0.0"));
		        	if (sampleState == i) buf.append(String.format("1.0"));
	            	buf.append(",");
		        }

	        	if (sampleState != states-1) buf.append(String.format("0.0"));
	        	if (sampleState == states-1) buf.append(String.format("1.0"));
		        
		        buf.append("}");
		        buf.append(",max" + type + "=");

		        buf.append(String.format("%d", sampleState ));
//		        buf.append(']'); 
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");

		        buf.append(String.format("%d", sampleState ));
		        buf.append(']');        		
        	}
        }
        
        buf.append(":");
        if (substitutions) {
            appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, node.getLength());
        }
        return buf.toString();
    }


    @Override
    public void close(PrintStream out) {
    	treeInput.get().close(out);
    }

	
}
