package beast.structuredCoalescent.math;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

import beast.core.Description;

/**
 * @author Nicola Felix Mueller
 */
@Description("Describes the derivates of the probabilities of configurations for ESCO as described in Mueller et al., 2016")
public class ode_integrator implements FirstOrderDifferentialEquations {

	double[] migration_rates;
	double[] coalescent_rates;
	double probs;
    int lineages;
    int states;
    int dimension;
    Integer[][] connectivity;
    Integer[][] sums;
    
    boolean belowzero = false;

    // constructor
    public ode_integrator(double[] migration_rates, double[] coalescent_rates, int lineages , int states, Integer[][] connectivity, Integer[][] sums){
        this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.lineages = lineages;
        this.states = states;
        belowzero = false;
        this.dimension = sums.length;
        this.connectivity = connectivity;
        this.sums = sums;
    }

    public int getDimension() {
        return dimension;
    }

    
    public void computeDerivatives(double t, double[] p, double[] pDot) {
    	// Calculates the change in probability of being in a state due to coalescence
    	for (int i = 0; i < p.length; i++){
    		pDot[i] = 0.0;
    		for (int s = 0; s < states; s++)
    			pDot[i] -= (sums[i][s]-1)*sums[i][s]*coalescent_rates[s];
    		
    		pDot[i] *= p[i];
    		// Stop the run if any configuration has a probability of lower than 0 of still existing
    		// normally the case when the integration time steps are too large
    		if (p[i]<0){
    			System.err.println("joint prob below 0");
    			System.exit(0);
    		}
    	}
    	// Calculate the change in the probability of being in a configuration due to migration
    	for (int i = 0; i < p.length; i++){
    		for (int j = 0; j < p.length; j++){
    			if (connectivity[i][j]!=null){
    				double m = p[i]*migration_rates[connectivity[i][j]];
    				pDot[i] -= m;
    				pDot[j] += m;
    			}
    		}
    	}
    }
        
    public static void main(String[] args) throws Exception{
        // 2d test
    	double[] migration_rates = {0.01, 0.001};
    	double[] coalescent_rates = {1.0, 1.0};
        int lineages = 2;
        int states = 2;
        /*
         * 0	1	1	0
         *	1	0	0	1
         *	1	0	0	1
         *	0	1	1	0
         */
        Integer[][] con = {{null,0,0,null},{1,null,null,0},{1,null,null,0},{null,1,1,null}};
        Integer[][] sums = {{2,0},{1,1},{1,1},{0,2}};

        FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.01);
        FirstOrderDifferentialEquations ode = new ode_integrator(migration_rates, coalescent_rates, lineages , states, con, sums);
        double[] y0 = new double[]{0,1,0,0};
        double[] y = new double[4];
    	integrator.integrate(ode, 0, y0, 100, y);

        System.out.println("Solution: " +y[0]+" "+y[1] + " " +y[2] + " " +y[3]);
    }

}