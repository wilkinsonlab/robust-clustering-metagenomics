# Robust Clustering Sub-states in Metagenomics

This pipeline, implemented in R, corresponds to an algorithm to automatically identify a reliable set of sub-states in longitudinal microbiome datasets. It is a generic and domain-independent procedure, applicable to whatever microbiome dataset.

The input are Operational Taxonomic Unit (OTU) vectors (in phyloseq or biom format, as it is explained below). The key consideration for grouping similar OTU vectors is according to Machine Learning clustering approaches, taking metagenomics beta diversity metrics as distance measures between samples. The procedure robustly defines the optimal number of clusters based on a comparison between several distance measures, distinct algorithms and different clustering scores. The variable factors that are combined to generate the robustness of our algorithm include: five different distance measures (JSD, rJSD, Bray-Curtis, Morisita-Horn and Kulczynski), two clustering scored (SI and Prediction Strength (PS) scores) followed by an additional bootstrapping process (evaluated with the Jaccard similarity score), and two distinct clustering approaches (PAM and Hclust). Additionally, different subsets of taxa could be considered ({all, dominant, non-dominant} x {species, genus level}.


## Basic Usage ##

Assuming all the required packages are already installed. You may need to download and install some missing dependencies.

From an interactive R session:
```r
# Load the functions
source('robust.clustering.metagenomics.functions.r')

# Example of call to the whole pipeline (with phyloseq object in a R data file):
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY")

# Or with .biom file:
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/David2014.biom",'David2014',"COLLECTION_DAY","~/RobustClustering/David2014/mapping_David2014.tsv")

# Example of call with dominant taxa at genus level (with phyloseq object in a R data file):
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY",'dominant','genus')

```

The file 'robust.clustering.metagenomics_callsExampleDatasets.r' includes cases of use of the different functions of the pipeline, both at high and low levels.


## Input and output description ##

  INPUTS:
  1. **path**: where the phyloseq object RData file is. 
  2. **RDataOrBiomFile**: the phyloseq object or the biom file.
    1. If it is a RData file, the phyloseq object must call 'data.norm' and it must have its OTU counts already normalized!!
    2. If it is a .biom file, the 5th argument is also required, because the .biom file will have only the OTU counts, which must be already normalized too!!
  3. **dataset.name**: used as label and suffix of files. Example: 'David2014', 'Chicks', etc.
  4. **variable.in.PCoAgraph**: name of a variable from sample_variables() in the phyloseq object, for the color of the samples in the PCoAgraph.
  5. [Optional: **taxaSubsetDominant**: string to determine if the taxa should be subsetted according to dominant taxa. Possible values: 'all' (default), 'dominant' or 'non-dominant'.]
  6. [Optional: **taxaSubsetGenus**: string to determine if taxa should be aggregated at genus level. Possible values: 'no' (default), 'genus' (compatible with dominant/non-dominant in taxaSubsetDominant parameter).]
  7. [Optional: **mapBiomFile**: mappping file instructions (only if biom format, else mapping info is included in the phyloseq object): comma separated values file, with samples in rows, being the first column the sampleID, and the remainder with their corresponding headers in the first row. The name to color the PCoA graphs must be one of these column headers within this mapping file.]

OUTPUTS:
  - phyloseq object with a new variable in the phyloseq object ($cluster) with the cluster identifier per sample. This object also is saved in 'data.normAndDist\_definitiveClustering\_\<dataset.label\>.RData'. It could be used as input of other R scripts with posterior steps of microbiome dynamics analysis.
  - The pipeline also returns several text and graph files with the results.



