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
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.IntervalType;
import beast.structuredCoalescent.math.independent_ode_integrator;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculate the probability of a tree under the structured coalescent assuming lineage independence with constant rates" +
		" as described in Mueller et. al.,2016. The Input rates are backwards in time migration rates" +
		" and pairwise coalescent rates translating to 1/Ne for m=1 in the Wright Fisher model")
public class IndependentStructuredCoalescent extends StructuredTreeDistribution implements Loggable {

    public Input<TraitSet> typeTraitInput = new Input<>(
            "typeTrait", "Type trait set.");     
    
    public Input<RealParameter> timeStepInput = new Input<RealParameter>(
    		"timeStep",
    		"how strong the seasonal patterns are, higher means less seasonal");

    public Input<RealParameter> coalescentRatesInput = new Input<RealParameter>(
    		"coalescentRates",
    		"how strong the seasonal patterns are, higher means less seasonal",Input.Validate.REQUIRED);  
    
    public Input<RealParameter> migrationRatesInput = new Input<RealParameter>(
    		"migrationRates",
    		"how strong the seasonal patterns are, higher means less seasonal",Input.Validate.REQUIRED);
    
    public Input<IntegerParameter> dim = new Input<IntegerParameter>(
    		"dim",
    		"how strong the seasonal patterns are, higher means less seasonal",Input.Validate.REQUIRED);  
    
    public Input<BooleanParameter> propTimeStepInput = new Input<BooleanParameter>(
    		"proportionalTimeStep",
    		"how strong the seasonal patterns are, higher means less seasonal"); 
        		
	
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] stateProbabilities;
	public DoubleMatrix[] stateFlow;
	
	private double max_posterior = Double.NEGATIVE_INFINITY;
	private double[] max_mig;
	private double[] max_coal;
    
    public int states;
    
    private boolean dependentHistory = true;
    private boolean traitInput = false;    
    private boolean normalize = true;
    private boolean propTimeStep = false;
    private int nr_lineages;    
    private double independent_likelihood=0.0;
    
    // Standard time Step of the RungeKutta if not specified different
    private double timeStep = 0.000001;
    
    
    // Set up for lineage state probabilities
    ArrayList<Integer> activeLineages;
    ArrayList<Double> lineStateProbs;
    
    @Override
    public void initAndValidate(){
    	// Check if the time steps should be done in such a way that they are proportional to the number
    	// of lineages present in a time interval
    	if (propTimeStepInput.get()!=null) propTimeStep = propTimeStepInput.get().getValue();
    	
    	// Calculate the tree intervals (time between events, which nodes participate at a event etc.)
    	treeIntervalsInput.get().calculateIntervals(); 	
       
        stateProbabilities = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];        
        stateFlow = new DoubleMatrix[2*treeIntervalsInput.get().getSampleCount() + 1];
        nrSamples = treeIntervalsInput.get().getSampleCount() + 1;
        
        double states_tmp = dim.get().getValue();
        states = (int) states_tmp;
        
        if (typeTraitInput.get() != null) traitInput = true;
        // Check all the boolean input parameters
        if (timeStepInput.get() != null) timeStep = timeStepInput.get().getValue();
        
        if(!dependentHistory) normalize = false;
        
        calculateLogP();
    }

    public double calculateLogP() {
    	// newly calculate tree intervals
    	treeIntervalsInput.get().calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	treeIntervalsInput.get().swap();
    	
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 
        lineStateProbs = new ArrayList<Double>();
        
        // Compute likelihood at each integration time and tree event starting at final sampling time and moving backwards
        logP = 0;

        // set the current time
        double currTime = 0.0;
        // total number of intervals
        final int intervalCount = treeIntervalsInput.get().getIntervalCount();
        // interval time counter
        int t = 0;
        // initialize the number of lineages
        nr_lineages = 0;
        // Captures the probabilities of lineages being in a state
        double[] p;	
        
        // Initialize the migration rates matrix
        double[][] migration_rates = new double[states][states];
        int c = 0;
        
       	for (int k = 0; k < states; k++) {
			for (int l = 0; l < states; l++){
				if (k!=l) {
					migration_rates[k][l] = migrationRatesInput.get().getArrayValue(c);
					c++;
				} 
				else{																// diagonal
					migration_rates[k][l] = 0.0;
				}
				
			}
		}
       	
       	// Initialize the coalescent rates
       	double[] coalescent_rates = new double[states];
		for (int k = 0; k < states; k++){
			coalescent_rates[k] = coalescentRatesInput.get().getArrayValue(k)/2;//(epiModelInput.get().getF(t,k,k) / (Y.get(k)*Y.get(k)));
		}

        
        // integrate until there are no more tree intervals
        do {
        	double nextIntervalTime = treeIntervalsInput.get().getInterval(t);

        	// Length of the current interval
        	final double duration = nextIntervalTime;// - currTime;
        	// if the current interval has a length greater than 0, integrate
        	if (duration > 0) {
        		if(dependentHistory)
        			p = new double[lineStateProbs.size()];		// Captures the probabilities of lineages being in a state
        		else
        			p = new double[lineStateProbs.size()+1];		// Captures the probabilities of lineages being in a state, last one keeps track of the probability
                
                // convert the array list to double[]
        		for (int i = 0; i<lineStateProbs.size(); i++)
                	p[i] = lineStateProbs.get(i);
                
        		// not needed
        		if(!dependentHistory)
        			p[lineStateProbs.size()] = 1;
        		
                double[] p_for_ode = new double[p.length];
                double ts = 0.0;
                
                // If proportial time step is true, set the integration time for the given interval 
                // inverse proportional to the number of lineages
                if (propTimeStep)
                	ts=timeStep/lineStateProbs.size();
                else 
                	ts=timeStep;
                
                // Never choose a longer time step than the integration window
                if(duration<(ts/2))
                	ts=duration/2;                
                
                FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(ts);
                // set the odes
                FirstOrderDifferentialEquations ode = new independent_ode_integrator(migration_rates, coalescent_rates, nr_lineages , states);
                // integrate	                
                integrator.integrate(ode, 0, p, duration, p_for_ode);
                
                // If the Dimension is larger than the maximum integer, at least one state prob is below 0 and the step is rejected
                if(ode.getDimension()==Integer.MAX_VALUE){
                	System.out.println(lineStateProbs.size());
                	System.out.println("lalalallal");
                	return Double.NEGATIVE_INFINITY;
                }
                
                for (int i = 0; i<lineStateProbs.size(); i++)
                	lineStateProbs.set(i, p_for_ode[i]);                
        	}
        	
        	// update the time
        	currTime = nextIntervalTime;
        	// event is coalescent event
        	if (treeIntervalsInput.get().getIntervalType(t) == IntervalType.COALESCENT){
        		logP += coalesce(t);
        		nr_lineages--;
	       	}       		

        	// event is sampling event
        	if (treeIntervalsInput.get().getIntervalType(t) == IntervalType.SAMPLE) {
       			addLineages(t);
       			nr_lineages++;
       		}    		
       		
        	// update the interval number
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
				for (int i = 0; i< states; i++)
					if (i == sampleState)
						lineStateProbs.add(1.0);
					else
						lineStateProbs.add(0.0);
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
				for (int i = 0; i< states; i++)
					if (i == sampleState)
						lineStateProbs.add(1.0);
					else
						lineStateProbs.add(0.0);
			}	
		}
    }
 
    
    private double coalesce(int currTreeInterval) {
    	double interval = 0.0;
    	
    	// normalize all lineage states in order to reduce the chance of too small numbers
    	// after multiplication of lineage states
    	if (normalize){    	
	    	for (int i = 0; i < nr_lineages; i++){
	    		double lineProbs = 0.0;
	    		for (int j = 0; j < states; j++)
	    			if (lineStateProbs.get(i*states+j)>=0.0){
	    				lineProbs += lineStateProbs.get(i*states+j);
	    			}else{
	    				System.out.println(lineStateProbs);
	    				System.err.println("state probability smaller than 0 before normalizing");
	    				System.exit(1);
	    			}
	    		for (int j = 0; j < states; j++)
	    			lineStateProbs.set(i*states+j, lineStateProbs.get(states*i+j)/lineProbs);   
	    		
	    		interval += Math.log(lineProbs);
	    	}
    	}
    	
    	// Get the daughter lineages
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
		
    	//  get the indices of the daughter lineages
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
		}
		DoubleMatrix lambda = DoubleMatrix.zeros(states);
		
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) { 
           	double coal = 0.0;
       		coal = coalescentRatesInput.get().getArrayValue(k)/2;
        	// calculate the coalscent rate. The 2 originates from lineage i coalescing with lineage j
       		// AND vice versa.
			final double pairCoalRate = coal * 2 *
					(lineStateProbs.get(daughterIndex1*states + k) * lineStateProbs.get(daughterIndex2*states + k));	
			
			if (!Double.isNaN(pairCoalRate)) lambda.put(k, pairCoalRate);
			else System.err.println("pairwise coalescent rate is nan");
        }
        
        activeLineages.add(coalLines.get(0).getParent().getNr());
                
		DoubleMatrix pVec = new DoubleMatrix();
		pVec.copy(lambda);
		if (normalize)
			pVec = pVec.div(pVec.sum());
				
		for (int i = 0; i < states; i++)
			lineStateProbs.add(pVec.get(i));

		stateProbabilities[coalLines.get(0).getParent().getNr() - nrSamples] = (pVec);		
        
		//Remove daughter lineages
		if (daughterIndex1>daughterIndex2){
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);
			for (int i = (states-1); i >= 0; i--)
				lineStateProbs.remove(daughterIndex1*states+i);
			for (int i = (states-1); i >= 0; i--)
				lineStateProbs.remove(daughterIndex2*states+i);
		}else{
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);
			for (int i = (states-1); i >= 0; i--)
				lineStateProbs.remove(daughterIndex2*states+i);			
			for (int i = (states-1); i >= 0; i--)
				lineStateProbs.remove(daughterIndex1*states+i);		
		}
		
     
		if (lambda.min()<0.0)
			System.err.println("Coalescent probability smaller than 0");
		
		if(!dependentHistory)
			interval = Math.log(independent_likelihood);
				
    	return Math.log(lambda.sum()) + interval;
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
     * Loggable interface
     */
    @Override
    public void init(PrintStream out)  {
    	for (int i=0;i<states;i++)
    		for (int j = 0; j<states;j++)
    			if(i!=j)
    				out.print("backward_rate_" + i +"to" +j+ "\t");
    	for (int i=0;i<states;i++)
    		for (int j = 0; j<states;j++)
    			if(i!=j)
    				out.print("forward_rate_" + j +"to" +i+ "\t");
    	for (int i=0;i<states;i++)
    		out.print("coalescent_rate_deme_" + i + "\t");
    	out.print("max_posterior\t");
    	for (int i=0;i<states*(states-1);i++)
    		out.print("max_mig_rate" + i + "\t");
    	for (int i=0;i<states;i++)
    		out.print("max_coal_rate" + i + "\t");

    }
    
    @Override
    public void log(int nSample, PrintStream out) {    	
    	int c = 0;
    	for (int i=0;i<states;i++)
			for (int j = 0; j<states;j++)
				if(i!=j){
					out.print(migrationRatesInput.get().getArrayValue(c) + "\t");
					c++;
				}

    	c = 0;
    	for (int i=0;i<states;i++)
    		for (int j = 0; j<states;j++)
    			if(i!=j){
    				out.print(migrationRatesInput.get().getArrayValue(c)/coalescentRatesInput.get().getArrayValue(j)*coalescentRatesInput.get().getArrayValue(i) + "\t");
    				c++;
    			}
    	for (int i=0;i<states;i++)
    		out.print(coalescentRatesInput.get().getArrayValue(i) +  "\t");
    	
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
