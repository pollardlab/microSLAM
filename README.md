<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/logo.png" width=200>

## Introduction

MicroSLAM is an R package to perform population structure leveraged association modeling for the microbiome from metagenomics data. Microbiome association studies typically link host disease or other traits to summary statistics measured in metagenomics data, such as diversity or taxonomic composition. But identifying disease-associated species based on their relative abundance does not provide insight into why these microbes act as disease markers, and it overlooks cases where disease risk is related to specific strains with unique biological functions. microSLAM is an implementation of a  mixed-effects model that performs association tests that connect host traits to the presence/absence of genes within each microbiome species, while accounting for strain genetic relatedness across hosts. Traits can be quantitative or binary (such as case/control).

MicroSLAM is fit in three steps for each species. The first step estimates a genetic relatedness matrix (GRM) that measures population structure of the microbial species across hosts. Step two calculates the association between population structure and the trait (y), enabling detection of species for which a subset of related strains confer risk. To identify specific genes (G) whose presence/absence across diverse strains is associated with the trait after adjusting for population structure, step three models the trait as a function of gene occurrence plus random effects (b) estimated from step two. Steps two and three can include adjustment for covariates measured on each host.

<p align="center">
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/Newflowchart.png" width=600>
</p>

### Below is a guide showing all operative functionality in R.
Install package in R

```
library(devtools)
install_github('pollardlab/microSLAM')
library(microSLAM)
```

### Inputs:

Inputs are usually two dataframes:

1) Gene data: a matrix containing information about the genes from a species' pangenome that are present versus absent in each sample (host), typically filtered to remove core genes (e.g., omiting genes that are more than 90% present across samples). Gene presence/absence can be estimated with bioinformatics tools, such as MIDAS. This will be a samples by genes matrix and must contain the column sample_name to indicate the sample name. The matrix provided in example_data is a simulated gene presence/absence matrix (.csv format) for 100 samples and 1000 genes.
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/genes.png" width=400>

2) Metadata: a matrix containing information about the phenotype (y) and the covariates for each sample. The sample names should match those in the gene data, and hence the metadata matrix may need to be pre-filtered on a species by species basis to contain data for only those samples with gene data for that species. The matrix provided in example_data (.csv format) is simulated metadata.
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/metadata.png" width=400>

Users will read in gene data from a file and metadata from a second file, both in .csv format, as shown below. The simulated example_data were modeled to have a strain that is correlated with the phenotype y, and a strain that is uncorrelated as well as three genes that are more associated with y than expected given the population structure. Age was randomly simulated as a covariate that is not associated with y.  

```
library(tidyverse)
library(magrittr)
data("exp_genedata") ### read in example gene data
## gene data is a number of samples by number of genes matrix for one species
data("exp_metadata") ### read in example metadata
## metadata is a number of samples (same as gene data) by k (number of covariates + 1) matrix
```

### Step 1: Calculate Genetic Relatedness Matrix (GRM)

Using the gene matrix imported above, step one will compute the GRM using the gene data. This GRM represents the population structure of the species across samples and is estimated as the similarity (1 minus Manhattan distance) of gene presence/absence vectors between pairs of samples. MicroSLAM will use the GRM to estimate sample-specific random effects for the mixed effects model. The gene matrix must contain the column sample_name for this step to run properly. Note that users may provide their own GRM in the following steps of microSLAM. For example, the GRM could be computed using a different distance or it could be estimated using other forms of genetic variation, such as core genome single-nucleotide variants.
```
GRM = calculate_grm(exp_genedata)
```
Example of GRM:

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/GRM.png" width=500>

#### Visualization of the GRM with the strain information used to generate it labeled
This gene data was generated with a strain (defined by a subset of correlated genes) that is associated with y in half of the samples and another strain that is not associated with y. In addition, three genes were simulated to be more related to y than is either strain.

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/exampleGRM.png" width=600>

For plotting GRM follow this code:

```
library(pheatmap)
myColor <- colorRampPalette(c("#FFFFFF","#009E73"))(50)
pheatmap(GRM,show_rownames=FALSE,show_colnames=FALSE,
    treeheight_row=0,treeheight_col = 0,
    labels_row="samples",labels_col="samples",
    main=paste("GRM"),
    border_color=NA,
    annotation_row = exp_metadata[,-5],
    color=myColor,
    clustering_distance_rows=as.dist(1-GRM),
    clustering_distance_col=as.dist(1-GRM),
    clustering_method="average")
```

### Step 2: Perform $\tau$ test for population structure (strain-trait associations)
Fit a baseline generalized linear model (glm) with only the covariate and an intercept, to obtain starting parameter estimates for the tau test. The family for the glm is binomial because y is binary.

```
glm_fit0 = glm("y ~ age + 1", data = exp_metadata, family = "binomial")
```

Fit a random effects glm using the baseline glm and GRM, and use it to estimate the parameter $\tau$, which measures how associated population structure is with the trait y. As long as there is a column sample_name in the metadata and the gene data, the program will check that the order is the same for GRM and the glm model, if not the user should know to check that the order of samples in the glm and the GRM should be the same.

```
glmm_fit = fit_tau_test(glm_fit0, GRM, species_id = "test", verbose = FALSE, log_file = NA)
summary.pop.struct.glmm(glmm_fit)
```

```
Species ID:  test
Formula:  y~age+1+b
family:  binomial logit
Fixed-effect covariates estimates:
 (Intercept) age
 -0.287 0.005
Converged:  TRUE
Number of iterations: 5
Tau:  2.354
Phi:  1 if logit or binomial should be 1
T value of tau: 0.624
Number of Samples: 100
```

This output can be interpreted as showing that microSLAM's glm in step two was able to converge, the $\tau$ variable is 2.354, the random effects variables (b) have been estimated, and the coefficients for the covariates have been estimated.

Next, test the significance of the estimated $\tau$ with a permutation test.

```
n_tau = 100 ### number of permutations for tau test
tautestfit = run_tau_test(glm_fit0, GRM, n_tau, species_id = "test", tau0=1, phi0=1)
```

Calculate the p-value from the permutation test.

```
pvalue = (sum(tautestfit$t >= glmm_fit$t) + 1)/n_tau
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/permutationnew.png" with=400>

In this case, the observed $\tau$ value (red vertical line) is larger than all $\tau$ values from 100 permutations that break the association between population structure and y (histogram shows permutation null distribution). This indicates that there is significant population structure for this species that is associated with y. To get a more precise p-value, more permutations could be run.  

Code to plot pvalues:

```
ggplot(tautestfit,aes(t))+
geom_histogram(bins = 30)+
geom_vline(xintercept = glmm_fit$t,color ="red")+
theme_minimal(base_size = 16)
```

### Step 3: $\beta$ test for gene-trait associations

Having detected population structure associated with the trait y, next fit a series of mixed effects models to test each gene for associations with y that are independent from the overall strain association. Rapidly gained and lost genes may have such independent associations. First, transform the gene data to a long matrix to fit each gene separately. If memory is a problem the gene data can be overwritten but we will use the square matrix to look into strain difference.

```
gene_long = exp_genedata %>%
    pivot_longer(cols = starts_with("gene"),
    names_to = "gene_id",
    values_to = "gene_value") %>%
    as.data.frame()
```

Fit a series of mixed effects glm models, one for each gene, that adjust for population structure using the random effects from Step two and perform a t-test to assess the significance of each gene's association with y ($\beta$).

```  
gene_test_df = fit_beta(glmm_fit,glm_fit0,GRM,gene_long,SPA=TRUE)
```

Plot a volcano plot of the results showing each gene's $\beta$ value versus its p-value. In this simulated data genes 1, 2, and 3 have been generated with associations that exceed the strain association. The genes with p-values of 0.005 or smaller are colored in red and correspond to these three genes.

```
ggplot(gene_test_df,aes(beta,-log10(SPA_pvalue)))+
 geom_point(data = gene_test_df[ which(gene_test_df$SPA_pvalue >= .005),], color = 'gray80', size=2)+
 geom_point(data = gene_test_df[which(gene_test_df$SPA_pvalue <= .005),], color = 'red', size = 2)+
 theme_minimal()
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/volcanoplotcolor.png?raw=true" with=400>


Example output dataframe for $\beta$ test

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/betadf1.png?raw=true">

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/betadf2.png?raw=true">


Columns are:

```
species_id: the id of the bacterial species
tau: estimate for tau genetic variance variable
gene_id: gene id for each gene
cor: correlation between the y variable and the gene data for that gene
cor_to_b: correlation between the random b variable and the gene data for that gene, which indicates genes driving overall strain-trait associations
z: z value estimated for gene from GLMM
var1: variance estimated for gene from GLMM
beta: beta estimate for gene from GLMM
se_beta: standard error for beta estimated from GLMM
t_adj: t value adjusted estimated from GLMM
SPA_pvalue: saddle point adjusted pvalue from GLMM
spa_score: saddle point adjusted t score from GLMM
SPA_zvalue:  saddle point adjusted z value from GLMM
pvalue_noadj: pvalue not adjusted from saddle point approximation
converged: did the saddle point adjustment (spa) algorithm converge or not
```

To visualize genes that contribute the most to the strain-trait association, each gene's correlation with the random effects vector b can be computed and the top correlated genes identified. In this example, genes with correlation above 0.4 are identified.

```
strain_genes = gene_test_df %>%
 filter(abs(cor_to_b) > .4)
rownames(exp_genedata) = exp_genedata$sample_name
strain_genes_values = exp_genedata[,strain_genes$gene_id]

pheatmap(strain_genes_values,
show_rownames=FALSE,
show_colnames=FALSE,
treeheight_row=0,
treeheight_col = 0,
labels_row = "samples",
labels_col = "samples",
main = paste("Strain Associated Genes"),
border_color = NA,
annotation_row = exp_metadata[,-5],
color = myColor,
clustering_method = "average")
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/strainheatmap.png?raw=true">

300 genes were modeled from the strain, and 288 of these have a correlation with b that exceeds 0.4. Samples could be clustered using these genes, and the genes can be considered as markers for the clades within the population structure of the species.
