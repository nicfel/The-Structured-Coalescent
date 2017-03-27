package beast.structuredCoalescent.math;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

import com.sun.org.apache.xalan.internal.xsltc.compiler.sym;

import beast.core.Description;


/**
 * @author Nicola Felix Mueller
 */
@Description("Describes the derivates of lineage state probabilities for LISCO as described in Mueller et al., 2016")
public class ode_masco implements FirstOrderDifferentialEquations {

	double[][] migration_rates;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
    
    boolean belowzero = false;

    // constructor
    public ode_masco(double[][] migration_rates, double[] coalescent_rates, int lineages , int states){
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        belowzero = false;
        this.dimension = this.lineages*this.states;
    }

    public int getDimension() {
        return dimension;
    }

    public void computeDerivatives(double t, double[] p, double[] pDot) {
    	// normalize each lineage such that the probability of the lineage not having coalesced
    	// is 1. This is done such that the coalescent rate is calculated conditional on all 
    	// lineages still existing/ not having coalesced
    	double[] p_norm = new double[p.length];
    	for (int i = 0; i<lineages; i++){
    		double linSum=0.0;
    		for (int j = 0; j<states; j++){
    			linSum += p[states*i+j];    			
    		}
    		for (int j = 0; j<states; j++){
    			p_norm[states*i+j] = p[states*i+j]/linSum;    
    		}
    	}
    	
    	// Compute the sum of line state probabilities for each state
    	int c = 0;
    	double[] sumStates = new double[states];
    	for (int j = 0; j<states; j++)
    		sumStates[j] = 0;
    	
    	for (int i = 0; i<lineages; i++){
    		for (int j = 0; j<states; j++){
    			sumStates[j] += p_norm[c];    			
    			c++;
    		}
    	}
    	
    	// Calculate dTdt
    	double dTdt = 0.0;
    	for (int i = 0; i<lineages; i++){
    		for (int j = 0; j<states; j++){
    			dTdt += p_norm[i*states+j]*coalescent_rates[j]*(sumStates[j] - p_norm[i*states+j]);
    		}
    	}
    	
    	// Caluclate the change in the lineage state probabilities for every lineage in every state
    	for (int i = 0; i<lineages; i++){
    		// get the change of P_t(T) due to coalescence w/o lineage i
    		double dTdtWoI = dTdt;
    		for (int j = 0; j<states; j++){
    			dTdtWoI -= 2* p_norm[i*states+j]*coalescent_rates[j]*(sumStates[j] - p_norm[i*states+j]);
    		}
//    		
//    		// calculate the rate at which coalescent events happen in the rest of the tree w/o lin i
//    		// calc sum state w/o i
//    		double[] sumStateWoi = new double[states];
//        	for (int a = 0; a<lineages; a++){
//        		if (a!=i){
//	        		for (int b = 0; b<states; b++){
//	        			sumStateWoi[b] += p_norm[a*states+b];    			
//	        		}
//        		}
//        	}        	
//        	
//    		double coalOut = 0.0;      	
//        	for (int a = 0; a<lineages; a++){
//        		if (a!=i){
//	        		for (int b = 0; b<states; b++){
//	        			coalOut += coalescent_rates[b] * p_norm[a*states+b]*(sumStateWoi[b]-p_norm[a*states+b]);    			
//	        		}
//        		}
//        	}   
//        	
//        	System.out.println(dTdtWoI + " " +  coalOut);
//        	System.out.println("");
        	        	
    		
    		// Calculate the probability of a lineage changing states
    		for (int j = 0; j<states; j++){
    			double migrates = 0.0;
    			for (int k = 0; k<states; k++){
    				if (j != k){
    					    					
    					// the probability of lineage i being in state j is p[i*nr_states +j]
    					migrates += p[states*i+k]*migration_rates[k][j] -
    									p[states*i+j]*migration_rates[j][k];
    					
    					// Check if the state probabilities do make sense, if the not,
    					// the integration time step is likely too large
    					if(p[states*i+k]<0.0){
    						// Use the dimension to communicate to the other side (the logP calculation)
    						dimension = Integer.MAX_VALUE;
    						return;
    					}   							
    				}
    			}// j    			 
    			
    			// Calculate the Derivate of p:
    			pDot[states*i+j] = migrates -
    					  2*coalescent_rates[j] * p[states*i+j] * (sumStates[j] - p_norm[states*i+j]) -
    					  	p[states*i+j]*dTdtWoI;
    		}// j
    	}// lineages    	
    }
    
    public static void main(String[] args) throws Exception{
        // 2d test
    	double[][] migration_rates = {{0.0, 0.001},{0.001, 0.0}};
    	double[] coalescent_rates = {0.01, 0.01};
        int lineages = 2;
        int states = 2;

        FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.01);
        FirstOrderDifferentialEquations ode = new independent_ode_integrator(migration_rates, coalescent_rates, lineages , states);
        double[] y0 = new double[]{1.,0.0,0.0,1.0};
        double[] y = new double[4];
    	integrator.integrate(ode, 0, y0, 100, y);
        System.out.println("Solution: " +y[0]+" "+y[1] + " " +y[2] + " " +y[3]);
        System.out.println("Solution Normalized: " +(y[0]/(y[0]+y[1]))+" "+(y[1]/(y[1]+y[0])) + " " +(y[2]/(y[2]+y[3])) + " "+  (y[3]/(y[2]+y[3])));
    }

}