## Simulations
This folder contains everything needed to recreate the results section in the text, as well as two additional folders not described in the main text. Some fodlers do not contain the exact output used in the paper, since the number of files resp. the size of those files is too large. They are however available upon request.

# Asymmetry


Contains all needed scripts to generate and analyze the results of the **Estimation of rate asymmetries** section. The number of files is too large to provide them on git.

# JointStateProbabilties

Contains matlab files to calculate and simulate the joint state probabilities for a two tip tree. This is in order to show how these joint state probabilities are the same when using the equtions described for ESCO and when directily simulating these joint state probabilities. This folder is not part of the paper.

# MigrationRates

Contains all needed scripts to generate and analyze the results of the **Estimation of migration rates** section. The relatively low number of files allows to put them of gitlab.

# OneState

This folder is not part of the paper, but it is an easy sanity check for the equations of MASCO, LISCO and SISCO as well as their implementation. It essentially contains files that compare the sampled tree heights under the standard unstructured coalescent model in BEAST2 as well as the sampled tree heights using MASCO, LISCO and SISCO for m=1

# RootStateProbabilities

Contains all needed scripts to generate and analyze the results of the **Root state probabilities** section. The relatively low number of files allows to put them of gitlab. The MultiTypeTree output tree file however is too large to be uploaded (over 100MB).

# SamplingBias

Contains all needed scripts to generate and analyze the results of the **Sampling bias** section. The number of files is too large (18000 files) to provide them on git.

# TreeHeightTest

Contains all needed scripts to generate and analyze the results of the **Tree height distributions under the structured co-alescent and its approximations** section. The relatively low number of files allows to put them of gitlab. Note that we use the standard BEAST2 operators to infer the tree heights using BASTA and not the MultiTypeTree operators. This we do, because the MultiTypeTree operators use the migration rates to propose new trees while we are interested in what the probability distribution is over trees and not how efficient the operators propose new steps

# Versions used

- BEAST 2.4.2 
- MASTER 5.0.2 
- BASTA 2.1.3 
- MuliTypeTree 6.2.1

# How to run the xmls

All the xmls with the ending \*esco.xml, \*masco.xml, \*lisco.xml and \*sisco.xml require the esco.jar from the /jar/ folder to run.

