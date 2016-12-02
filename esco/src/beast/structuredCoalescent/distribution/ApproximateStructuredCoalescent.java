/*
 * Copyright (C) 2016 Nicola Felix Mueller (nicola.felix.mueller@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package beast.structuredCoalescent.distribution;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.IntervalType;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculate thes probability of a tree under the structured coalescent with constant rates" +
			" as described in Volz et. al.,2012. The Input rates are backwards in time migration rates" +
			" and pairwise coalescent rates translating to 1/Ne for m=1 in the Wright Fisher model")
public class ApproximateStructuredCoalescent extends StructuredTreeDistribution implements Loggable {

	
    public Input<RealParameter> coalescentRatesInput = new Input<RealParameter>(
    		"coalescentRates",
    		"how strong the seasonal patterns are, higher means less seasonal"); 
    
    public Input<RealParameter> migrationRatesInput = new Input<RealParameter>(
    		"migrationRates",
    		"how strong the seasonal patterns are, higher means less seasonal");
    
    public Input<TraitSet> typeTraitInput = new Input<>(
            "typeTrait", "Type trait set.");    
    
    public Input<RealParameter> timeStepInput = new Input<RealParameter>(
    		"timeStep",
    		"how strong the seasonal patterns are, higher means less seasonal");
    
    public Input<IntegerParameter> dim = new Input<IntegerParameter>(
    		"dim",
    		"how strong the seasonal patterns are, higher means less seasonal",Input.Validate.REQUIRED);  

	
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] stateProbabilities;
	private DoubleMatrix preMigration;
    
    public int states;
    
    private boolean traitInput = false;
    private boolean preMigrationFilled = false;
        
    // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
    ArrayList<Integer> activeLineages;
    ArrayList<DoubleMatrix> lineStateProbs;
    
    @Override
    public void initAndValidate() {
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)
    	treeIntervalsInput.get().calculateIntervals();
 	        
        stateProbabilities = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];        
        nrSamples = treeIntervalsInput.get().getSampleCount() + 1;

        // direct conversion to integer didn't seem to work, so take rout via double
        double states_tmp = dim.get().getValue();
        states = (int) states_tmp;        
        
        // check if there is a set of trait values given as Input (untested)
        if (typeTraitInput.get() != null) traitInput = true;
        
        // Calculate the marginal likelihood
        calculateLogP();
    }

    public double calculateLogP() {
    	
    	// Initialize preMigration. preMigration saves the result of the matrix exponentiation
    	// for the timeStep given in the input. This saves time to recalculate the matrix
    	// exponential for the input time step
    	preMigration = new DoubleMatrix(states,states);
    	preMigrationFilled = false;
    	
    	// newly calculate tree intervals
    	treeIntervalsInput.get().calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	treeIntervalsInput.get().swap();
    	        
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 
        lineStateProbs = new ArrayList<DoubleMatrix>();
        
        // Compute likelihood at each integration time and tree event starting at final sampling time and moving backwards
        ArrayList<Double> integrationTimes = new ArrayList<Double>();
        double time = 0.0;
        while (time <= treeIntervalsInput.get().getTotalDuration()){
        	integrationTimes.add(time);
        	time += timeStepInput.get().getValue();
        }
        integrationTimes.add(time);
        time += timeStepInput.get().getValue();
        integrationTimes.add(time);
        
        // what tree interval are we in?
        int currTreeInterval = 0;
        // current time (since present)
        double currTime = integrationTimes.get(0);
        // time of next integration step
        double nextIntegrationTime;
        // time next tree interval starts
        double nextIntervalTime = treeIntervalsInput.get().getInterval(0);
        final int intervalCount = treeIntervalsInput.get().getIntervalCount();
        logP = 0;
        // interval time counter
        int t = 0;
                
        
        do {       
        	nextIntegrationTime = integrationTimes.get(t+1);													// Update Integration time 
        	while (nextIntervalTime <= nextIntegrationTime) {													// while there are still events to take place before the next integration step
        		/*
        		 * compute contribution of last interval to likelihood
        		 */
	        	final double duration = nextIntervalTime - currTime;
	        	if (duration > 0) {
				  	
	        		// Update line states
	        		boolean negInf = updateLineProbs(t, duration); 
	  				if (negInf){
	  					logP = Double.NEGATIVE_INFINITY;	        		
	  				}
	        		// Calculate interval contribution
	        		logP += -computeLambdaSum(t,duration) * duration;
	        	}                
	               
	        	/*
	        	 *  compute contribution of event to likelihood
	        	 */
	        	if (treeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.COALESCENT) {	//CHECK IF INDEX IS CORRECT HERE
	        		final double pairCoalRate = computeLambdaPair(t, currTreeInterval);
	        		if (pairCoalRate > 0) {
	        			logP += Math.log(pairCoalRate);
	        		}else{
	        			System.out.println("probability of a coalescent event is smaller than zero");
	        			return Double.NEGATIVE_INFINITY;
	        		}      		
	        		updateLineProbsCoalEvent(t, currTreeInterval);												// Set parent lineage state probs and remove children
		       	}
	       		
	        	/*
	        	 * add new lineage
	        	 */
	       		if (treeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.SAMPLE) {
	       			addLineages(currTreeInterval);
	       		}
	       		
	       		currTime += duration;
	       		currTreeInterval++; 
	       		
	       		if (currTreeInterval < intervalCount) {
	       			nextIntervalTime = currTime + treeIntervalsInput.get().getInterval(currTreeInterval);
        		} else {
	       			nextIntervalTime = Double.POSITIVE_INFINITY;
	       			currTreeInterval--; //stay in last interval
	       		}

        	}
        	final double duration = nextIntegrationTime - currTime;
	    	if (duration > 0) { 
        		// Update line states
        		boolean negInf = updateLineProbs(t, duration); 
  				if (negInf){
  					logP = Double.NEGATIVE_INFINITY;	        		
  				}
        		// Calculate interval contribution
        		double add = -computeLambdaSum(t,duration) * duration;
        		logP += add;
    		}		    	
    		currTime += duration;
    		t++;
    		
        }while(integrationTimes.get(t-1) <= treeIntervalsInput.get().getTotalDuration());
	             
        //Compute likelihood of remaining tree intervals (coal events occuring before origin)
        if (Double.isInfinite(logP))logP = Double.NEGATIVE_INFINITY;
        return logP;
   	
    }

        
    
    /*
     * Should be correct
     */
    private double computeLambdaPair(int t, int currTreeInterval) {
		
		List<Node> coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
    	if (coalLines.size() > 2) {
			System.err.println("Unsupported coalescent at non-binary node");
			System.exit(0);
		}
    	if (coalLines.size() < 2) {
    		System.out.println();
    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
    		System.out.println();
    		return Double.NaN;
		}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			return 0.0;
		}
		double lambda = 0.0;
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) { 
        	final double popCoalRate = coalescentRatesInput.get().getArrayValue(k);
        	// divide by 2 in order to have the pairwise coalescent rates as input
			final double pairCoalRate = popCoalRate * 
				2*lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(k)/2;
			
			if (!Double.isNaN(pairCoalRate)) lambda += pairCoalRate;
			
    	}

    	return lambda;
    }
    
    /*
     * Should be correct, matrix stuff not checked
     */
    private double computeLambdaSum(int t, double dt) {
    	
		int lineCount = activeLineages.size();
		double lambdaSum = 0.0; 										// holds the sum of the pairwise coalescent rates over all lineage pairs
		
		if (lineCount < 2) return lambdaSum;											// Zero probability of two lineages coalescing	


		/*
		 * Sum over line pairs (scales better with numbers of demes and lineages)
		 */	
		DoubleMatrix popCoalRates;
		/*
		 * If the transmission matrix has only elements in the diagonal
		 * do fancy matrix calculation epiModelInput.get().diagTrans
		 */
		popCoalRates = new DoubleMatrix(states);
		for (int k = 0; k < states; k++){
        	// divide by 2 in order to have the pairwise coalescent rates as input
			popCoalRates.put(k, coalescentRatesInput.get().getArrayValue(k)/2);
		}
		
		DoubleMatrix sumStates = DoubleMatrix.zeros(states);
		DoubleMatrix diagElements = DoubleMatrix.zeros(states);
		for (int linI = 0; linI < lineCount; linI++) {
			sumStates = sumStates.add(lineStateProbs.get(linI));
			diagElements = diagElements.add(lineStateProbs.get(linI).mul(lineStateProbs.get(linI)));
		}
		DoubleMatrix M = sumStates.mul(sumStates);
		lambdaSum = M.sub(diagElements).mul(popCoalRates).sum();
			
    	return lambdaSum;    	
    }
    
    private boolean updateLineProbs(int t, double dt) {    	

		/*
		 * Calculate the rate at which states change. State changes in time t
		 * happen from k to l. The total rate at which a state changes depends on
		 * the migration (G-matrix) and the transmission form deme k to l (F-matrix)
		 */
    	DoubleMatrix Q = DoubleMatrix.zeros(states,states);
    	
    	int c = 0;
   	
		for (int k = 0; k < states; k++) 
		{
			for (int l = 0; l < states; l++) 
			{
				if (k != l) 
				{																        				// off-diagonal
					final double totalRate = migrationRatesInput.get().getArrayValue(c);
					c++;
					Q.put(k, l, totalRate);
					Q.put(k, k, Q.get(k,k)-totalRate);
				}
			}
		}

		/*
		 * Update lineage state probabilities for all active lineages. The change in state
		 * probability is calculated as: i) Q the state change matrix per time is multiplied 
		 * with time ii) Q is multiplied with the current state probabilities such that
		 * dP = Q*dt*s has the elements [q11 * s1 , q12 * s1...; q21 *s2...].
		 * The changes in state probabilities is than the sum of columns.  
		 */		
			
		if (dt == timeStepInput.get().getValue() && !preMigrationFilled)
		{
			preMigration = matrixExponential(Q, dt);
			DoubleMatrix newP = new DoubleMatrix();
			newP.copy(preMigration);
			for (int lin = 0; lin < activeLineages.size(); lin++) {
				DoubleMatrix oldProbs = lineStateProbs.get(lin);			
				DoubleMatrix probs = vecMat(newP, oldProbs);			
				/*
				 * divide the state probabilities by the sum of state probabilities
				 * in order to avoid errors due to numerical errors as small as they
				 * might be
				 */
				probs = probs.div(probs.sum());
				lineStateProbs.set(lin, probs);
			}
			preMigrationFilled=true;
		}
		else 
			{
			if (dt == timeStepInput.get().getValue() && preMigrationFilled)
			{
				DoubleMatrix newP = new DoubleMatrix();
				newP.copy(preMigration);
				for (int lin = 0; lin < activeLineages.size(); lin++) {
					DoubleMatrix oldProbs = lineStateProbs.get(lin);			
					DoubleMatrix probs = vecMat(newP, oldProbs);			
					/*
					 * divide the state probabilities by the sum of state probabilities
					 * in order to avoid errors due to numerical errors as small as they
					 * might be
					 */
					probs = probs.div(probs.sum());
					lineStateProbs.set(lin, probs);
				}
			}
			else
			{
				DoubleMatrix newP = matrixExponential(Q, dt);
				for (int lin = 0; lin < activeLineages.size(); lin++) {
					DoubleMatrix oldProbs = lineStateProbs.get(lin);			
					DoubleMatrix probs = vecMat(newP, oldProbs);			
					/*
					 * divide the state probabilities by the sum of state probabilities
					 * in order to avoid errors due to numerical errors as small as they
					 * might be
					 */
					probs = probs.div(probs.sum());
					lineStateProbs.set(lin, probs);
				}
			}
		}			


		
		/*
		 * Everything seems to have worked correctly
		 */
		return false;
		
    }
    
   

    
    private void updateLineProbsCoalEvent(int t, int currTreeInterval) {
    	List<Node> coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
//    	if (coalLines.size() > 2) {
//			throw new Exception("Unsupported coalescent at non-binary node");
//		}
		int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		
		if (daughterIndex1 == -1 || daughterIndex2 == -1){
	    	treeIntervalsInput.get().swap();
	    	coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
	    	daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
	    	daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
//	    	if (daughterIndex1 == -1 || daughterIndex2 == -1)
//				throw new Exception("Active lineages does not contain coalescing lineages");
		}
		List<Node> parentLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		
		
		//Add parent to activeLineage and initialize parent's state prob vector
		Node parentNode = parentLines.get(0);
		activeLineages.add(parentNode.getNr());
		DoubleMatrix pVec = new DoubleMatrix(states);

		//Compute parent lineage state probabilities in pVec
		DoubleMatrix coalRates = DoubleMatrix.zeros(states, states);
		double lambda;
		for (int k = 0; k < states; k++) {
    		lambda = coalescentRatesInput.get().getArrayValue(k) 
    			*2*lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(k);
        	
        	coalRates.put(k, k, lambda);            
        }	
		/*
		 * pVec := new state probabilities of the lineage dt after coalescing		
		 */
		pVec = coalRates.rowSums().div(coalRates.sum());
        lineStateProbs.add(pVec);
		stateProbabilities[parentNode.getNr() - nrSamples] = (pVec);
        
		//Remove daughter lineages
		activeLineages.remove(daughterIndex1);
		lineStateProbs.remove(daughterIndex1);
		
		daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr()); //index may have changed after removal of first daughter
		activeLineages.remove(daughterIndex2);
		lineStateProbs.remove(daughterIndex2);
    }
    
    
    
    private void addLineages(int currTreeInterval) {
		List<Node> incomingLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		if(traitInput){
			/*
			 * If there is a typeTrait given as Input the model will take this
			 * trait as states for the taxons
			 */		
			for (Node l : incomingLines) {				
				activeLineages.add(l.getNr());
				int sampleState = (int) typeTraitInput.get().getValue(l.getID());
				DoubleMatrix sVec = DoubleMatrix.zeros(states);
				sVec.put(sampleState, 1.0);			
				lineStateProbs.add(sVec);
			}			
		}else{		
			/*
			 * If there is no trait given as Input, the model will simply assume that
			 * the last value of the taxon name, the last value after a _, is an integer
			 * that gives the type of that taxon
			 */
			for (Node l : incomingLines) {
				activeLineages.add(l.getNr());
				String sampleID = l.getID();
				String[] splits = sampleID.split("_");
				int sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
				DoubleMatrix sVec = DoubleMatrix.zeros(states);
				sVec.put(sampleState, 1.0);			
				lineStateProbs.add(sVec);
			}	
		}
    }
        
    public DoubleMatrix getStateProb(int nr){
    	return stateProbabilities[nr - nrSamples];
    }
    
    public DoubleMatrix[] getStateProbabilities(){
    	return stateProbabilities;
    }
    
    public String getType(){
    	if (typeTraitInput.get()!=null){
    		return typeTraitInput.get().getTraitName();
    	}
    	else{
    		return "state";
    	}
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	return true;
    }
    
    /*
     * Quickly added vector times matrix 
     */
    private DoubleMatrix vecMat(DoubleMatrix mat, DoubleMatrix vec){
    	DoubleMatrix prod = DoubleMatrix.zeros(states);
    	for (int i = 0; i < states; i++)
    		for (int j = 0; j < states; j++)
    			prod.put(i, prod.get(i) + vec.get(j)*mat.get(j,i));
    	return prod;
    }
    
    private DoubleMatrix matrixExponential(DoubleMatrix Q, double dt){
    	return MatrixFunctions.expm(Q.mul(dt));
    }
    
    
  	/*
     * Loggable interface
     */
//    @Override
//    public void init(PrintStream out)  {
//		out.print("Event\tInterval\t");
//    }
//    
//    @Override
//    public void log(int nSample, PrintStream out) {
//		out.print(event + "\t");
//		out.print(interval + "\t");
//    }
//    
//    @Override
//    public void close(PrintStream out) {
//    }
    
    
}