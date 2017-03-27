### How to prepare your data
All the structured coalescent methods in the package esco are using the sequence names to determine where lineages were sampled. In order for the methods to know where the lineages were sampled, their names have to be changed slightly. Add an underscore + a number indicating the sampling state, for example:

~~~
KJ189303|1998|Colombia -> KJ189303|1998|Colombia_0
GQ868561|1999|Colombia -> GQ868561|1999|Colombia_0
JF804021|1998|PuertoRico -> JF804021|1998|PuertoRico_1
...
~~~

Note that java starts indexing at 0 and **not** 1. The sequences will now be in the location or state determined by the number after the underscore.

### How to prepare your xml
An xml that can be used with any of the methods can be created easily starting from an xml with a **constant coalescent** prior (or any other population prior for that matter). Create an xml in BEAUTi with an alignment that has the sequence names as described above. After the sampling time, site model, clock rates, the coalescent prior etc. are set up, save the xml. To get the xml running follow the next few steps:

**First step:** open the newly created xml and replace the first line with:

```xml
<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate='Standard' beautistatus='' namespace="beast.core:beast.evolution.alignment:beast.structuredCoalescent.distribution:beast.structuredCoalescent.logger:beast.structuredCoalescent.operator:beast.structuredCoalescent.model:beast.structuredCoalescent.operators:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" version="2.0">
```

This needs to be done in order to tell beast which java classes it needs to load

**Second step:** Now the initial values for the backwards migration rates and for the coalescent rates need to be defined. This needs to be done such that the parameters are inferable. Look for the following line in the xml:

```xml
</state>
```

Above that line, add the migration and coalescent rates by copying the lines below just above the line </state>

```xml
<!-- Lisco parameters -->
<parameter id="coalRates" name="stateNode">0.2 0.2 0.1</parameter>
<parameter id="migRates" name="stateNode">0.01 0.02 0.001 0.02 0.1 0.01</parameter>
```

The coalRates are the pairwise coalescent rates, i.e. the rate at which two lineages in the same state coalesce. The number of entries (here 0.2 0.2 0.1) need to be equal to the number of different states. migRates are the backwards migration rate. Their number of entries (here 0.01 0.02 0.001 0.02 0.1 0.01) needs to be equal to nr states*(nr states - 1)

**Third step:** Now, the xml needs to be changed such that the structure coalescent is used as a population prior. To do so, look for 4 lines in the xml (they should be a couple of lines below the parameters specified before) that look something like:

```xml
<distribution id="CoalescentConstant.t:denv1_truncated_1" spec="Coalescent">
    <populationModel id="ConstantPopulation.t:denv1_truncated_1" spec="ConstantPopulation" popSize="@popSize.t:denv1_truncated_1"/>
    <treeIntervals id="TreeIntervals.t:denv1_truncated_1" spec="TreeIntervals" tree="@Tree.t:denv1"/>
</distribution>
```

the id names can slightly deviate from the example here, but if you look for *<distribution id="CoalescentConstant.t*, you should always find the correct line, in case the initial xml was setup with a constant coalescent prior. These lines tell Beast to use an unstructured constant coalescent prior. Replace those lines with:

```xml
<distribution id="structuredCoalescent" spec="util.CompoundDistribution">
    <distribution id="coalescent" spec="ExactStructuredCoalescent" dim="3">
        <structuredTreeIntervals spec="StructuredTreeIntervals" id="TreeIntervals" tree="@treename"/>
        <coalescentRates idref="coalRates"/>
        <migrationRates idref="migRates"/>
        <parameter id="stepSize" name="timeStep">4</parameter>
    </distribution>
</distribution>
<distribution spec='beast.math.distributions.Prior' x="@coalRates">
    <distr spec="beast.math.distributions.Exponential" mean="1"/>
</distribution>
<distribution spec='beast.math.distributions.Prior' x="@migRates">
    <distr spec="beast.math.distributions.Exponential" mean="1"/>
</distribution>
```

It is important that the section where it says tree="@treename" has the same treename for the newly copied lines as it had for the constant coalescent. In this example we would need to change the things as follow

~~~
 tree="@treename" -> tree="@Tree.t:denv1"
~~~

the entry *spec="IndependentStructuredCoalescent"* tells BEAST to use ESCO as a coalescent prior. To switch between the structured coalescent and its approximations, replace *ExactStructuredCoalescent* with:

~~~
ExactStructuredCoalescent (for ESCO)
Masco (for MASCO)
IndependentStructuredCoalescent (for LISCO)
ApproximateStructuredCoalescent (for SISCO)
~~~

Next, check that the number of dimensions is correct it should say dim="nr states"

The lines containing ...Exponential... are prior distributions on the coalescent rates and migration rates with a mean=1 here, which can be changed.

**Fourth step:** Next, operators to scale the migration and coalescent rates need to be added. Replace the line that has the word PopSizeScaler in it, e.g.:

```xml
<operator id="PopSizeScaler.t:denv1_truncated_1,2" spec="ScaleOperator" parameter="@popSize.t:denv1_truncated_1,2" scaleFactor="0.75" weight="3.0"/>
```

with

```xml
<operator id="CoalescentScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="false" parameter="@coalRates" weight="1.0"/>

<operator id="MigrationScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="false" parameter="@migRates" weight="1.0"/>
```

This allows the migration and coalescent rates to be scaled by BEAST. Replacing the PopSizeScale is not actually neccessary, but elsewise, beast would some times update the effective population size which does not influence the posterior in any way, since we replace the constant coalescent population prior

**Fifth step:** Now the only thing left is to log the parameters and the trees. To log the migration and coalescent rates, look for a line that looks something like:

```xml
<logger id="tracelog" fileName="denv1_lisco.log" logEvery="50000" model="@posterior" sanitiseHeaders="true" sort="smart">
```

add the following below that line:

```xml
<!-- Lisco-->
<log idref="migRates"/>
<log idref="coalRates"/>
<log idref="structuredCoalescent"/>
```

This now logs the migration rates, the coalescent rates and the structured coalescent probability. The only thing left to do is to log the trees. Look for a couple of lines that look like:

```xml
<logger id="treelog.t:denv1_truncated_1,2" fileName="denv1_lisco.trees" logEvery="50000" mode="tree">
    <log id="TreeWithMetaDataLogger.t:denv1_truncated_1,2" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:denv1"/>
</logger>
```

They will be almost at the end of the file. Add the following lines below them

```xml
<logger id="probtreelog" fileName="denv1_lisco.states.trees" logEvery="50000" mode="tree">
    <log id="logTrees" spec="StructuredTreeLogger" tree="@treename" exactDensity="@coalescent"/>
</logger>
```

Again, it is important that the section where it says tree="@treename" has the same treename for the newly copied lines as it had for the constant coalescent. In this example we would need to change the things as follows:

~~~
tree="@treename" -> tree="@Tree.t:denv1"
~~~

Depending on which coalescent approach is defined as a structured coalescent with the id="coalescent", the input of the tree logger needs to be adjusted and *exactDensity* needs to be changed as follows

~~~
ESCO -> exactDensity
MASCO -> mascoDensity
LISCO -> independentDensity
SISCO -> approximateDensity
~~~

The trees that are logged this way are using the maximum state probabilities **conditioned only on the subtree below every node** as colors. They also log the probability of every node being in any state conditioned on the subtree but not the full tree.
