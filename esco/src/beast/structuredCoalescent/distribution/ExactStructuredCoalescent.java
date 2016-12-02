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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.jblas.DoubleMatrix;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.IntervalType;
import beast.structuredCoalescent.math.ode_integrator;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculate the probability of a tree under the exact numerical structured coalescent with constant rates" +
		" as described in Mueller et. al.,2016. The Input rates are backwards in time migration rates" +
		" and pairwise coalescent rates translating to 1/Ne for m=1 in the Wright Fisher model")
public class ExactStructuredCoalescent extends StructuredTreeDistribution implements Loggable {
	
    public Input<TraitSet> typeTraitInput = new Input<>(
            "typeTrait", "Type trait set.");        
  
    public Input<RealParameter> timeStepInput = new Input<RealParameter>(
    		"timeStep",
    		"the time step used for rk4 integration");

    public Input<RealParameter> coalescentRatesInput = new Input<RealParameter>(
    		"coalescentRates",
    		"pairwise coalescent rates that translate to 1/Ne in the Wright Fisher model",Input.Validate.REQUIRED); 
    
    public Input<RealParameter> migrationRatesInput = new Input<RealParameter>(
    		"migrationRates",
    		"Backwards in time migration rate",Input.Validate.REQUIRED);
    
    public Input<IntegerParameter> dim = new Input<IntegerParameter>(
    		"dim",
    		"the number of different states",Input.Validate.REQUIRED);  
        		
	
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] nodeStateProbabilities;
	public ArrayList<Double> jointStateProbabilities;
	public ArrayList<Integer> numberOfLineages;
	public Integer[][] connectivity;
	public Integer[][] sums;
	public Integer[] sumsTot;
	public ArrayList<ArrayList<Integer>> combination;
	
	private double[] migration_rates;
	private int[][] migration_map;
	private double[] coalescent_rates;
			
	
	private double max_posterior = Double.NEGATIVE_INFINITY;
	private double[] max_mig;
	private double[] max_coal;
    
    public int states;
    
    private boolean traitInput = false;    
    private int nr_lineages;
        
    // Standard time Step of the RungeKutta if not specified different
    private double timeStep = 0.000001;
    
    
    // Set up for lineage state probabilities
    ArrayList<Integer> activeLineages;
    ArrayList<Double> lineStateProbs;
    
    @Override
    public void initAndValidate(){
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)
    	treeIntervalsInput.get().calculateIntervals(); 	
       
        nodeStateProbabilities = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];        
        nrSamples = treeIntervalsInput.get().getSampleCount() + 1;
        
        // direct conversion to integer didn't seem to work, so take rout via double
        double states_tmp = dim.get().getValue();
        states = (int) states_tmp;        
        
        // check if there is a set of trait values given as Input (untested)
        if (typeTraitInput.get() != null) traitInput = true;
        if (timeStepInput.get() != null) timeStep = timeStepInput.get().getValue();
        
        migration_rates = new double[states*(states-1)];
        migration_map = new int[states][states];
        coalescent_rates = new double[states];
        
        // Calculate the marginal likelihood
        calculateLogP();
    }

    public double calculateLogP() {
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)    	
    	treeIntervalsInput.get().calculateIntervals();
    	treeIntervalsInput.get().swap();
    	
        // Set up for lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 
        lineStateProbs = new ArrayList<Double>();
        
        // Compute likelihood at each integration time and tree event starting at final sampling time and moving backwards
        logP = 0;          
        
        // Initialize the line state probabilities
        // total number of intervals
        final int intervalCount = treeIntervalsInput.get().getIntervalCount();
        // counts in which interval we are in
        int t = 0;
        nr_lineages = 0;
        // Captures the probabilities of lineages being in a state
        double[] p;			
        
        // Initialize the migration rates matrix
        int c = 0;
        
       	for (int k = 0; k < states; k++) {
			for (int l = 0; l < states; l++){
				if (k!=l) {
					migration_rates[c] = migrationRatesInput.get().getArrayValue(c);
					migration_map[k][l] = c;
					c++;
				}
				else{
					coalescent_rates[k] = coalescentRatesInput.get().getArrayValue(k)/2;
				}
				
			}
		}

        boolean first = true;
        // integrate until there are no more tree intervals
        do {
        	double nextIntervalTime = treeIntervalsInput.get().getInterval(t);
        	// Length of the current interval
        	final double duration = nextIntervalTime;// - currTime;
        	// if the current interval has a length greater than 0, integrate
        	if (duration > 0) {
    			p = new double[jointStateProbabilities.size()];		// Captures the probabilities of lineages being in a state
                
                // convert the array list to double[]
        		for (int i = 0; i<jointStateProbabilities.size(); i++)
                	p[i] = jointStateProbabilities.get(i);                
        		
                double[] p_for_ode = new double[p.length];
                
                double ts=timeStep;
                if(duration<timeStep)
                	ts=duration/2;

                // initialize integrator
                FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(ts);
                // set the odes
                FirstOrderDifferentialEquations ode = new ode_integrator(migration_rates, coalescent_rates, nr_lineages , states, connectivity, sums);
                // integrate	                
                integrator.integrate(ode, 0, p, duration, p_for_ode);
                
                // if the dimension is equal to the max integer, this means that a calculation
                // of a probability of a configuration resulted in a value below 0 and the
                // run will be stopped
                if(ode.getDimension()==Integer.MAX_VALUE){
                	System.out.println("lalalallal");
                	return Double.NEGATIVE_INFINITY;
                }
                
                // set the probabilities of the system being in a configuration again
                for (int i = 0; i<p_for_ode.length; i++)
                	jointStateProbabilities.set(i, p_for_ode[i]);
                
        	}
        	/*
        	 *  compute contribution of event to likelihood
        	 */
        	if (treeIntervalsInput.get().getIntervalType(t) == IntervalType.COALESCENT){
        		nr_lineages--;
        		logP += coalesce(t);
	       	}
       		
        	/*
        	 * add new lineage
        	 */
       		if (treeIntervalsInput.get().getIntervalType(t) == IntervalType.SAMPLE) {
       			nr_lineages++;
       			addLineages(t, first);
       			first = false;
       		}
       		
       		
       		t++;
        }while(t < intervalCount);
	             
        //Compute likelihood of remaining tree intervals (coal events occuring before origin)
        if (Double.isInfinite(logP))logP = Double.NEGATIVE_INFINITY;
        if (max_posterior<logP && logP<0){
        	max_posterior=logP;
        	max_mig = new double[states*(states-1)];
        	max_coal = new double[states];
        	for (int i = 0 ; i < 1;i++)
        		max_mig[i] = migrationRatesInput.get().getArrayValue(i);
        	for (int i = 0; i < 1; i++)
        		max_coal[i] = coalescentRatesInput.get().getArrayValue(i);
        }
     
        return logP;   	
    }

        
    
    private void addLineages(int currTreeInterval, boolean first) {
		List<Node> incomingLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		int sampleState=0;
		if(traitInput){
			/*
			 * If there is a typeTrait given as Input the model will take this
			 * trait as states for the taxons
			 */		
			for (Node l : incomingLines) {				
				activeLineages.add(l.getNr());
				sampleState = (int) typeTraitInput.get().getValue(l.getID());
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
				sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
			}	
		}
		
		int nrs = 0;
		if (first)
			nrs = states;
		else
			nrs = combination.size()*states;
		
		Integer[][] newSums = new Integer[nrs][states];
		Integer[] newSumsTot = new Integer[nrs];
		
		
		// Add all new combinations of lineages
		ArrayList<Double> newJointStateProbabilities = new ArrayList<Double>();
		ArrayList<ArrayList<Integer>> newCombination = new ArrayList<ArrayList<Integer>>();
		if (first){
			for (int i = 0; i < states; i++){
				ArrayList<Integer> add = new ArrayList<Integer>();
				add.add(i);
				if (i == sampleState)
					newJointStateProbabilities.add(1.0);
				else
					newJointStateProbabilities.add(0.0);
				newCombination.add(add);
			}
			
			// get the state in which lineages (n=1 in this case) are
			for (int i = 0; i < states; i++){
				newSumsTot[i] = 0;
				for (int j = 0; j < states; j++){
					if (j == i)
						newSums[i][j] = 1;
					else
						newSums[i][j] = 0;

				}
			}
			
		}
		else {
			// Dublicate all entries by a factor of #states
			for (int i = 0; i < combination.size(); i++){
				ArrayList<Integer> addBase = combination.get(i);
				for (int j = 0; j < states; j++){						
					ArrayList<Integer> newComb = new ArrayList<Integer>(addBase);
					newComb.add(j);
					newCombination.add(newComb);
					if (j == sampleState)
						newJointStateProbabilities.add(jointStateProbabilities.get(i));
					else
						newJointStateProbabilities.add(0.0);
				}
				// calculate the new number of lineages in a state for each configuration
				for (int j = 0; j < states; j++){
					newSums[states*i+j][j] = sums[i][j]+1;
					for (int l = 0; l < states; l++){
						if (l != j)
							newSums[states*i+j][l] = sums[i][l];
					}
				}
			}
			
			// Get the number of pairs of lineages in each configuration
			for (int i = 0; i < nrs; i++){
				int news = 0;				
				for (int j = 0; j < states; j++){
					double add = newSums[i][j]-1;
					if (add>0)
						news += add;
				}
				newSumsTot[i] = news;
			}
			
		}
		// set the old values (pre event) to the new ones
		combination = newCombination;	
		jointStateProbabilities = newJointStateProbabilities;
		
		// newly build the connectivity matrix (how to transition between configurations
		connectivity = new Integer[combination.size()][combination.size()];
		// build the connectivity matrix
		for (int a = 0; a < combination.size(); a++){
			for (int b = 0; b < combination.size(); b++){
				int diff = 0;
				int[] directs = new int[2];
				ArrayList<Integer> comb1 = combination.get(a);
				ArrayList<Integer> comb2 = combination.get(b);

				for (int i = 0; i < comb1.size(); i++){
					int d = comb1.get(i) - comb2.get(i);
					if (d != 0){
						diff++;
						directs[0] = comb1.get(i);
						directs[1] = comb2.get(i);
					}
				}
				if (diff == 1){
					connectivity[a][b] = migration_map[directs[0]][directs[1]];
				}
			}
		}
		sums = newSums;
		sumsTot = newSumsTot;
		
    }
 
    
    private double coalesce(int currTreeInterval) {
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
    	
    	// get the indices of the two daughter lineages
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
		}
		
		// check which index is large such that the removing starts 
		// with the one with the larger value
		if (daughterIndex1>daughterIndex2){
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);
		}else{
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);
		}

		// add the new parent lineage as an active lineage
		activeLineages.add(coalLines.get(0).getParent().getNr());
				
		// calculate the number of combinations after the coalescent event
		int nrs = combination.size()/states;
		
		// newly initialize the number of lineages per configuration
		Integer[][] newSums = new Integer[nrs][states];
		Integer[] newSumsTot = new Integer[nrs];

		
		// find all joint probabilities where the two lineages are in the same deme
		ArrayList<Double> newProbability = new ArrayList<Double>();
		ArrayList<ArrayList<Integer>> newCombination = new ArrayList<ArrayList<Integer>>();
		double[] pairwiseCoalRate = new double[states];
		int futureState = 0;
		for (int i = 0; i < jointStateProbabilities.size(); i++){
			// Checks if it is a configuration where both daughter lineages are in the same state
			if (combination.get(i).get(daughterIndex1) == combination.get(i).get(daughterIndex2)){
				ArrayList<Integer> coalLoc = new ArrayList<Integer>(combination.get(i));
				
				newSums[futureState] = sums[i];	
				newSums[futureState][combination.get(i).get(daughterIndex1)]--;
				futureState++;
				
				if (daughterIndex1>daughterIndex2){
					coalLoc.remove(daughterIndex1);
					coalLoc.remove(daughterIndex2);
				}else{
					coalLoc.remove(daughterIndex2);
					coalLoc.remove(daughterIndex1);
				}
				coalLoc.add(combination.get(i).get(daughterIndex1));
				newCombination.add(coalLoc);	
				
				newProbability.add(coalescent_rates[combination.get(i).get(daughterIndex1)]*jointStateProbabilities.get(i));
				pairwiseCoalRate[combination.get(i).get(daughterIndex1)] += 2*coalescent_rates[combination.get(i).get(daughterIndex1)]*jointStateProbabilities.get(i);
			}
		}
		
		combination = newCombination;
		jointStateProbabilities = newProbability;
		
		connectivity = new Integer[combination.size()][combination.size()];
		// build the connectivity matrix
		for (int a = 0; a < combination.size(); a++){
			for (int b = 0; b < combination.size(); b++){
				int diff = 0;
				int[] directs = new int[2];
				ArrayList<Integer> comb1 = combination.get(a);
				ArrayList<Integer> comb2 = combination.get(b);

				for (int i = 0; i < comb1.size(); i++){
					int d = comb1.get(i) - comb2.get(i);
					if (d != 0){
						diff++;
						directs[0] = comb1.get(i);
						directs[1] = comb2.get(i);
					}
				}
				if (diff == 1){
					connectivity[a][b] = migration_map[directs[0]][directs[1]];
				}
			}
		}
		
		
		
		for (int i = 0; i < nrs; i++){
			int news = 0;				
			for (int j = 0; j < states; j++){
				int add = newSums[i][j]-1;
				if (add>0)
					news += add;
			}
			newSumsTot[i] = news;
		}

		

	
		// do normalization
		double prob = 0.0;
		for (int i = 0; i < pairwiseCoalRate.length; i++)
			prob += pairwiseCoalRate[i];
		
		for (int i = 0; i < jointStateProbabilities.size(); i++)
			jointStateProbabilities.set(i,jointStateProbabilities.get(i)/prob);
		
		DoubleMatrix pVec = new DoubleMatrix(states);
		
		
		for (int i = 0; i < pairwiseCoalRate.length; i++)
			pVec.put(i, pairwiseCoalRate[i]/prob);
		
		nodeStateProbabilities[coalLines.get(0).getParent().getNr() - nrSamples] = pVec;
		
		sums = newSums;
		sumsTot = newSumsTot;
		
		// return the normlization constant as a probability (in log space)
    	return Math.log(prob);
    }
    
    public DoubleMatrix getStateProb(int nr){
    	return nodeStateProbabilities[nr - nrSamples];
    }
    public DoubleMatrix[] getStateProbabilities(){
    	return nodeStateProbabilities;
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
     * Loggable interface
     */
    @Override
    public void init(PrintStream out)  {
    	out.print("max_posterior\t");
    	for (int i=0;i<states*(states-1);i++)
    		out.print("max_mig_rate" + i + "\t");
    	for (int i=0;i<states;i++)
    		out.print("max_coal_rate" + i + "\t");

    }
    
    @Override
    public void log(int nSample, PrintStream out) {    	
    	out.print(max_posterior +"\t");
    	for (int i=0;i<states*(states-1);i++)
    		out.print(max_mig[i] + "\t");
    	for (int i=0;i<states;i++)
    		out.print(max_coal[i] + "\t");

    }
    
    @Override
    public void close(PrintStream out) {
    }
    
}
