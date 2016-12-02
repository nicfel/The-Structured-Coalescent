/*
* File MCMC.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package beast.structuredCoalescent.operator;

import java.io.IOException;
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.HashSet;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
//import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.util.Randomizer;

@Description("MCMC chain. This is the main element that controls which posterior " +
        "to calculate, how long to run the chain and all other properties, " +
        "which operators to apply on the state space and where to log results.")
@Citation(value=
        "Bouckaert RR, Heled J, Kuehnert D, Vaughan TG, Wu C-H, Xie D, Suchard MA,\n" +
                "  Rambaut A, Drummond AJ (2014) BEAST 2: A software platform for Bayesian\n" +
                "  evolutionary analysis. PLoS Computational Biology 10(4): e1003537"
        , year = 2014, firstAuthorSurname = "bouckaert",
        DOI="10.1371/journal.pcbi.1003537")
// basically the original MCMC class that allows for fixing the order of the operators. This is needed to have some sort
// of step wise increasement of rate while always newly inferring the migration histories for x-times.
// Is only used for inferring the root states using MultiTypeTree
public class MCMC_fixedOperators extends MCMC {
	
    public Input<BooleanParameter> alwaysAcceptInput = new Input<>(
    		"alwaysAccept", "if true, a new proposed step will always be accepted (Default false)");
    
    public Input<BooleanParameter> operatorsFixedInput = new Input<>(
    		"operatorsFixed", "if true,the operator schedule will always be the same (Default false)");

    
//    private boolean alwaysAccept = false;
    private boolean operatorsFixed = true;
    
    // list with the weights of the operators
    private int[] operatorWeights;
    private int[] operatorWeightsRemaining;
    private int nrOperators;
//    private int currOperator = 0;
    
    
    
    public MCMC_fixedOperators() {
    	System.out.println("lsdldfsaldfsaldafsljkagdskjhfgskhjgalhjkadghjklghjkdghjkdg");
    }


    @Override
    public void initAndValidate() {
//    	if(alwaysAcceptInput.get()!=null) alwaysAccept = alwaysAcceptInput.get().getValue();
//    	if(operatorsFixedInput.get()!=null) operatorsFixed = alwaysAcceptInput.get().getValue();
    	
//    	System.exit(0);
    	nrOperators = operatorsInput.get().size();
    	
    	if(operatorsFixed){
    		operatorWeights = new int[operatorsInput.get().size()];
    		operatorWeightsRemaining = new int[operatorsInput.get().size()];
    		for (int i = 0; i < operatorsInput.get().size(); i++){
    			operatorWeights[i] = (int) operatorsInput.get().get(i).getWeight();
    			operatorWeightsRemaining[i] = (int) operatorsInput.get().get(i).getWeight();
    			System.out.println(operatorWeights[i]);
    		}
//    		System.exit(0);
    	}

    	super.initAndValidate();
    } // init

    public void log(final int sampleNr) {
    	super.log(sampleNr);
    } // log

    public void close() {
        super.close();
    } // close

//    protected double logAlpha;
//    protected boolean debugFlag;
//    protected double oldLogLikelihood;
//    protected double newLogLikelihood;
//    protected int burnIn;
//    protected int chainLength;
//    protected Distribution posterior;

    protected List<Logger> loggers;

    @Override
    public void run() throws IOException, SAXException, ParserConfigurationException {
    	super.run();
    } // run;

    @Override
    protected void doLoop() throws IOException {
        int corrections = 0;
        final boolean isStochastic = posterior.isStochastic();
        
        if (burnIn > 0) {
        	Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
        }
        for (int sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {
            final int currentState = sampleNr;

            state.store(currentState);
//            if (m_nStoreEvery > 0 && sample % m_nStoreEvery == 0 && sample > 0) {
//                state.storeToFile(sample);
//            	operatorSchedule.storeToFile();
//            }
            Operator operator = operatorsInput.get().get(0);
            for (int i = 0; i<nrOperators;i++ ){
            	if(operatorWeightsRemaining[i]>0){
            		operator = operatorsInput.get().get(i);
            		operatorWeightsRemaining[i]--;
            		break;
            	}
            	if(operatorWeightsRemaining[i]==0 && i==(nrOperators-1)){
            		for (int j = 0; j<nrOperators;j++ ){
            			operatorWeightsRemaining[j] = operatorWeights[j];
            		}
            		operator = operatorsInput.get().get(0);
            		operatorWeightsRemaining[0]--;
//            		for (int j = 0; j<nrOperators;j++ ){
//            			System.out.println(operatorWeightsRemaining);
//            		}
//                    System.exit(0);
            	}
            }
            	
//    		for (int j = 0; j<nrOperators;j++ ){
//    			System.out.println(operatorWeightsRemaining);
//    		}
//
//            
//            System.out.println(operator);
//            System.exit(0);


            final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
            Evaluator evaluator = null;

            if (evaluatorDistribution != null) {
                evaluator = new Evaluator() {
                    @Override
                    public double evaluate() {
                        double logP = 0.0;

                        state.storeCalculationNodes();
                        state.checkCalculationNodesDirtiness();

                        try {
                            logP = evaluatorDistribution.calculateLogP();
                        } catch (Exception e) {
                            e.printStackTrace();
                            System.exit(1);
                        }

                        state.restore();
                        state.store(currentState);

                        return logP;
                    }
                };
            }
            final double logHastingsRatio = operator.proposal(evaluator);

            if (logHastingsRatio != Double.NEGATIVE_INFINITY) {

            	if (operator.requiresStateInitialisation()) {
            		state.storeCalculationNodes();
            		state.checkCalculationNodesDirtiness();
            	}

                newLogLikelihood = posterior.calculateLogP();

                logAlpha = newLogLikelihood - oldLogLikelihood + logHastingsRatio; //CHECK HASTINGS

                if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
                    // accept
                    oldLogLikelihood = newLogLikelihood;
                    state.acceptCalculationNodes();

                    if (sampleNr >= 0) {
                        operator.accept();
                    }
                } else {
                    // reject
                    if (sampleNr >= 0) {
                        operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
                    }
                    state.restore();
                    state.restoreCalculationNodes();
                }
                state.setEverythingDirty(false);
            } else {
                // operation failed
                if (sampleNr >= 0) {
                    operator.reject(-2);
                }
                state.restore();
				if (!operator.requiresStateInitialisation()) {
                    state.setEverythingDirty(false);
                    state.restoreCalculationNodes();
				}
            }
            log(sampleNr);

//            if (debugFlag && sampleNr % 3 == 0 || sampleNr % 10000 == 0) {
//                // check that the posterior is correctly calculated at every third
//                // sample, as long as we are in debug mode
//            	final double originalLogP = isStochastic ? posterior.getNonStochasticLogP() : oldLogLikelihood;
//                final double logLikelihood = isStochastic ? state.robustlyCalcNonStochasticPosterior(posterior) : state.robustlyCalcPosterior(posterior);
//                if (isTooDifferent(logLikelihood, originalLogP)) {
//                    reportLogLikelihoods(posterior, "");
//                    Log.err.println("At sample " + sampleNr + "\nLikelihood incorrectly calculated: " + originalLogP + " != " + logLikelihood
//                    		+ "(" + (originalLogP - logLikelihood) + ")"
//                            + " Operator: " + operator.getClass().getName());
//                }
//                if (sampleNr > NR_OF_DEBUG_SAMPLES * 3) {
//                    // switch off debug mode once a sufficient large sample is checked
//                    debugFlag = false;
//                    if (isTooDifferent(logLikelihood, originalLogP)) {
//                        // incorrect calculation outside debug period.
//                        // This happens infrequently enough that it should repair itself after a robust posterior calculation
//                        corrections++;
//                        if (corrections > 100) {
//                            // after 100 repairs, there must be something seriously wrong with the implementation
//                        	Log.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
//                            state.storeToFile(sampleNr);
//                            operatorSchedule.storeToFile();
//                            System.exit(1);
//                        }
//                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);;
//                    }
//                } else {
//                    if (isTooDifferent(logLikelihood, originalLogP)) {
//                        // halt due to incorrect posterior during intial debug period
//                        state.storeToFile(sampleNr);
//                        operatorSchedule.storeToFile();
//                        System.exit(1);
//                    }
//                }
//            } else {
//                if (sampleNr >= 0) {
//                	operator.optimize(logAlpha);
//                }
//            }
            callUserFunction(sampleNr);

            // make sure we always save just before exiting
            if (storeEvery > 0 && (sampleNr + 1) % storeEvery == 0 || sampleNr == chainLength) {
                /*final double logLikelihood = */
                state.robustlyCalcNonStochasticPosterior(posterior);
                state.storeToFile(sampleNr);
                operatorSchedule.storeToFile();
            }
        }
        if (corrections > 0) {
        	Log.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
        }
    }
    
} // class MCMC
