# The Structured Coalescent and its Approximations

Nicola F. Müller<sup>1,2</sup>, David A. Rasmussen<sup>1,2</sup>, Tanja Stadler<sup>1,2</sup>

<sup>1</sup>ETH Zurich, Department of Biosystems Science and Engineering, 4058 Basel, Switzerland

<sup>2</sup>Swiss Institute of Bioinformatics (SIB), Switzerland




## Abstract
Phylogenetics can be used to elucidate the movement of pathogens between different host populations when the location of samples are considered alongside of pathogen sequence data. Pathogen phylogenies therefore offer insights into the movement of pathogens not available from classic epidemiological data alone.
However, current phylogeographic methods to quantify migration patterns from phylogenies have several known shortcomings. 
In particular, one of the most widely used method treats migration the same as mutation, and as such does not incorporate information about population demography.
This may lead to severe biases in estimated migration rates for datasets where sampling is biased across populations. 
On the other hand, the structured coalescent allows us to coherently model the migration and transmission process, but current implementations struggle with complex datasets due to the need to additionally infer ancestral migration histories.
Thus, approximations to the structured coalescent which integrate over all ancestral migration histories have been developed. However, the validity and robustness of these approximations remain unclear.
We here provide an exact numerical solution to the structured coalescent that does not require the inference of migration histories. 
While this solution is computationally unfeasible for large datasets, it clarifies the assumptions of previously developed approximate methods and allows us to provide an improved approximation to the structured coalescent.
We have implemented these methods in BEAST2, and we show how our newly described approach outperforms previously described methods in accuracy at comparable computational cost.

## License

The content of this project itself is licensed under the [Creative Commons Attribution 3.0 license](http://creativecommons.org/licenses/by/3.0/us/deed.en_US), and the java source code of **esco** is licensed under the [GNU General Public License](http://www.gnu.org/licenses/).