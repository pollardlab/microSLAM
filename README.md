<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/logo.png" width=200>

## Introduction

microSLAM is an R package to perform population structure leveraged association modeling for the microbiome from metagenomics data. 
Microbiome association studies typically link host disease or other traits to summary statistics measured in metagenomics data, 
such as diversity or taxonomic composition. 
But identifying disease-associated species based on their relative abundance does not provide insight into why these microbes act as disease markers, 
and it overlooks cases where disease risk is related to specific strains with unique biological functions. 
microSLAM is an implementation of a  mixed-effects model that performs association tests that connect host traits to the 
presence/absence of genes within each species, while accounting for strain genetic relatedness across hosts. 
Traits can be quantitative or binary (such as case/control).

microSLAM is fit in three steps for each species. 
The first step estimates a genetic relatedness matrix (GRM) that measures population structure of the microbial species across hosts. 
Step two calculates the association between population structure and the trait (y), enabling detection of species for which a subset of related strains confer risk. To identify specific genes (G) whose presence/absence across diverse strains is associated with the trait after adjusting for population structure, step three models the trait as a function of gene occurrence plus random effects (b) estimated from step two. Steps two and three can include adjustment for covariates measured on each host.

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

### Input:

Inputs are usually two data frames per species:

1) Sample-by-Gene Binary Data: a binary matrix representing the presence or absence of genes from a species' pangenome across multiple samples.
This matrix is typically filtered to remove core genes, such as those present in more than 90% of samples. 
Gene presence/absence can be assessed using metagenotyping tools like [MIDAS](https://github.com/czbiohub-sf/MIDAS). 
The data must follow a samples-by-genes format and must include a `sample_name` as the first column to represent sample identifiers.
The example dataset provided in example_data is a simulated gene presence/absence matrix in CSV format, containing 100 samples and 1,000 genes.
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/genes.png" width=400>

2) Sample Metadata: a matrix containing phenotype data (`y`) and other covariates for each sample.
The sample names must match those in the gene data, ensuring consistency between datasets.
Since the list of samples may vary across species, the metadata must be filtered on a species-by-species basis to ensure it includes only samples that have corresponding gene data.
The example dataset provided in the example_data is a simulated metadata matrix in CSV format.
<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/metadata.png" width=400>


Users will import gene data and metadata from separate CSV files, as demonstrated below.

```
library(tidyverse)
library(magrittr)
data("exp_genedata")
data("exp_metadata")
```

### Step 1: Calculate Genetic Relatedness Matrix (GRM)

Using the imported sample-by-gene matrix, the first step is to compute the genetic relatedness matrix (GRM) based on the gene presence/absence data.
This GRM captures the population structure of the species across samples and 
is calculated as the pairwise similarity score, defined as 1 minus the Manhattan distance between gene presence/absence vectors. 
microSLAM utilized the GRM to estimate sample-specific random effects within the mixed effects model. 
To ensure proper execution, the gene matrix must include a sample_name column.
Alternatively, users can provide their own GRM, e.g. using a different distance metric.

```
GRM = calculate_grm(exp_genedata)
```

Example of GRM:

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/GRM.png" width=500>

#### Visualization of the GRM with the strain information used to generate it labeled

The simulated example_data were designed to include:

- A strain that is correlated with the phenotype `y`
- A strain that is uncorrelated with `y`
- Three genes exhibiting a stronger association with `y` than expected give population structure.
- Age, a randomly simulated as a covariate that is not associated with `y`.  

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
    annotation_row = exp_metadata[,2:5],
    color=myColor,
    clustering_distance_rows=as.dist(1-GRM),
    clustering_distance_col=as.dist(1-GRM),
    clustering_method="average")
```

### Step 2: $\tau$ test for strain-trait association

Fit a baseline generalized linear model (GLM) that includes only the covariate and an intercept to establish 
initial parameter estimates for the tau test. Since `y` is binary, the GLM uses a binomial family.

```
glm_fit0 = glm("y ~ age + 1", data = exp_metadata, family = "binomial")
```

Fit a random-effects generalized linear model (GLM) using the baseline GLM and the Genetic Relatedness Matrix (GRM) to estimate the parameter $\tau$, 
which quantifies the association between population structure and the trait (`y`).
To ensure proper alignment, verify that the sample order in the GRM matches that in the GLM model before procedding.

```
glmm_fit=fit_tau_test(glm_fit0, GRM, species_id="species_test", verbose = FALSE, log_file="./microSLAM.log")
summary.pop.struct.glmm(glmm_fit)
```

```
Species ID:  test
Formula:  y~age+1+b
family:  binomial logit
Fixed-effect covariates estimates:
 (Intercept) age
 -0.306 0.005
Converged:  TRUE
Number of iterations: 5
Tau:  2.525
Phi:  1 if logit or binomial should be 1
T value of tau: 0.667
Number of Samples: 100
```

This output indicates that microSLAM's GLM in step two successfully converged, with the $\tau$ parameter estimated at 2.525. 
Additionally, the random effects variables (`b`) have been computed, and the coefficients for the covariates have been estimated.

Next, test the significance of the estimated $\tau$ with a permutation test.

```
n_tau = 100 
tautestfit = run_tau_test(glm_fit0, GRM, n_tau, species_id = "species_test", tau0=1, phi0=1, seed=63)
```

Calculate the p-value from the permutation test.

```
pvalue = (sum(tautestfit$t >= glmm_fit$t) + 1)/n_tau
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/permutationnew.png" with=400>

In this case, the observed $\tau$ value (represented by the red vertical line) is greater than 
all $\tau$ values derived from 100 permutations, where the association between population and `y` was randomized. 
The histogram illustrated that the null distribution of $\tau$ values from these permutations.
This result suggests that population structure is significantly associated with y for this species. 
To obtain a more precise p-value, additional permutations can be conducted.

R code to plot p-values from the permutation test:

```
ggplot(tautestfit,aes(t))+
  geom_histogram(bins = 30)+
  geom_vline(xintercept = glmm_fit$t,color ="red")+
  theme_minimal(base_size = 16)
```

### Step 3: $\beta$ test for gene-trait associations

After detecting population structure associated with the trait (`y`), 
the next step is to fit a series of mixed effects models to evaluate independent associations
between each gene and `y`,
while adjusting for population structure using the random effects estimated in Step 2.
A t-testis then performed to assess the significance of each gene's association with y ($\beta$).

Genes that are rapidly gained or lost may exhibit such independent associations.

```  
gene_test_df = fit_beta(glmm_fit, glm_fit0, GRM, exp_genedata, SPA=TRUE)
```

Generate a volcano plot displaying each gene’s $\beta$ value against its p-value. 
In this simulated dataset, genes 1, 2, and 3 were designed to exhibit associations that surpass the strain association. 
Genes with p-values ≤ 0.005 are highlighted in red, corresponding to these three significant genes.

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
species_id: Identifier of the bacterial species
tau: Estimate for tau genetic variance variable
gene_id: Unique identifier for each gene.
cor: Correlation between the response variable (y) and the gene presence/absence data
cor_to_b: Correlation between the random effects variable (b) and the gene data, indicating genes that contribute tp overall strain-trait associations
z: Z-score estimate for each gene from the GLMM
var1: Variance estimate for each gene from the GLMM
beta: Regressuin coefficient (beta) estimate for each gene from the GLMM
se_beta: Standard error for beta estimate from the GLMM
t_adj: Adjusted t value estimate from the GLMM
SPA_pvalue: Saddlepoint-adjusted pvalue from the GLMM
spa_score: Saddlepoint-adjusted t-score from the GLMM
SPA_zvalue: Saddle point adjusted z-score from the GLMM
pvalue_noadj: Unadjusted pvalue
converged: Indicates whether the SPA algorithm sucessfully converted or not
```

To identify genes that contribute most to the strain-trait association, 
the correlation between each gene and the random effects vector (b) was computed. 
The top correlated genes can then be selected based on their correlation values. 
In this example, genes with a correlation greater than 0.4 are plotted.

```
strain_genes = gene_test_df %>%
 filter(abs(cor_to_b) > .4)
  rownames(exp_genedata) = exp_genedata$sample_name
  strain_genes_values = exp_genedata[,strain_genes$gene_id]

pheatmap(strain_genes_values,
  show_rownames=FALSE, show_colnames=FALSE,
  treeheight_row=0, treeheight_col = 0,
  labels_row = "samples", labels_col = "samples",
  main = paste("Strain Associated Genes"),
  border_color = NA,
  annotation_row = exp_metadata[,2:5],
  color = myColor,
  clustering_method = "average")
```

<img src="https://github.com/miriam-goldman/microSLAM/blob/main/other/strainheatmap.png?raw=true">

A total of 300 genes were modeled from the strain, with 288 genes exhibiting a correlation greater than 0.4 with the random effects vector (b). 
These genes could be used to cluster samples, serving as markers for clades within the population structure of the species.

## Run the full microSLAM pipeline

We also offer a wrapper function that allows users to execute all the above steps in a single workflow, 
including GRM computation, GLM fitting, tau estimation, permutation testing, and gene association testing.

We strongly recommend that users customize the `run_microslam` function to suit their specific analysis needs.

```
run_microslam(exp_genedata, exp_metadata, data.frame(gene_id = colnames(exp_genedata)[-1]), 
              GRM, "species_test", "./example_output", "./example_log", do_SPA = TRUE, 
              formula_string = "y ~ age + 1", response_var = "y", 
              family = "binomial", n_tau = 100, verbose = TRUE)
```

