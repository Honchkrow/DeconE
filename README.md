- [Section 1: Introduction](#section-1-introduction)
- [Section 2: Installation](#section-2-installation)
- [Section 3: Gnerating Pseudo Bulk Data](#section-3-gnerating-pseudo-bulk-data)
  - [Section 3.1: Gnerating Pseudo Bulk Data From Massive RNA-seq Studies](#section-31-gnerating-pseudo-bulk-data-from-massive-rna-seq-studies)
  - [Section 3.2: Gnerating Pseudo Bulk Data From scRNA-seq data](#section-32-gnerating-pseudo-bulk-data-from-scrna-seq-data)
  - [Code Demo 1: Evaluating The Deconvolution Results In A Simple Manner](#code-demo-1-evaluating-the-deconvolution-results-in-a-simple-manner)
- [Section 4: Noise Analysis](#section-4-noise-analysis)
  - [Section 4.1: Generating Noised Bulk Data With Different Models](#section-41-generating-noised-bulk-data-with-different-models)
  - [Code Demo 2: Evaluating The Deconvolution Results For Multiple Method](#code-demo-2-evaluating-the-deconvolution-results-for-multiple-method)
- [Section 5: Rare Component Analysis](#section-5-rare-component-analysis)
- [Section 6: Single Cell Related Functions](#section-6-single-cell-related-functions)
- [Section 7: Well-Characterized Deconvolution Datasets](#section-7-well-characterized-deconvolution-datasets)
- [Citation](#citation)

**Links**:

- [Deconer manual](https://honchkrow.github.io/Deconer/inst/documents/Deconer_manual.pdf)
- [Deconer vignettes](https://honchkrow.github.io/Deconer/inst/documents/Deconer_intro.html)
- [Deconer datasets](https://honchkrow.github.io/Deconer_dataset/)

## Section 1: Introduction

Cell type proportion is related to phenotypes or diseases ([Wang, *et al.*](https://doi.org/10.1038/s41467-018-08023-x), [Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)). Therefore, quantifying cell or tissue proportions is important for understanding the mechanisms in biological processes.

Here, we propose a cell type deconvolution evaluating toolkit named '**Deconer**' to perform comprehensive and systematic evaluation of different algorithms.

**Deconer** consists of 6 main functional components as described below:

- Pseudo bulk data generation (from large-scale bulk data and single-cell data).
- Stability analysis under various types of noise.
- Rare component analysis.
- Unknown component analysis.
- Comprehensive evaluation metrics, along with visually appealing figure generation.
- Well-characterized datasets for deconvolution utilities.

In the subsequent sections, we will provide a comprehensive tutorial on how to use Deconer.

## Section 2: Installation

Deconer is based on R and can be easily installed on Windows, Linux as well as MAC OS.

First, users should install [R >= 4.1.0](https://www.r-project.org/).

Next, install devtools and Deconer.

```R
# install devtools
install.packages('devtools')

# install the Deconer package
devtools::install_github('Honchkrow/Deconer')

# load Deconer
library(Deconer)
```

## Section 3: Gnerating Pseudo Bulk Data

Generating pseudo-bulk data is a challenging problem. Deconer takes inspiration from previous works ([Francisco, *et al.*](https://doi.org/10.1038/s41467-020-20288-9), [Wang, *et al.*](https://doi.org/10.1038/s41467-018-08023-x), [Wei, *et al.*](https://doi.org/10.1093/bib/bbab362), [Racle, *et al.*](https://doi.org/10.7554/eLife.26476)) and the [Tumor Deconvolution DREAM Challenge](https://www.synapse.org/#!Synapse:syn15589870/wiki/) to provide various strategies for generating pseudo-bulk RNA-seq data from both large-scale bulk RNA-seq data and scRNA-seq data.

### Section 3.1: Gnerating Pseudo Bulk Data From Massive RNA-seq Studies

Many deconvolution methods require cell type-specific bulk data as prior knowledge during the deconvolution process. To create a more realistic simulation, we have collected 302 well-characterized bulk RNA-seq data for generating pseudo-bulk data. During the data generation process, one-third of the samples are used for generating the mixture, while the remaining samples are used for generating the external reference.

Taking inspiration from the [Tumor Deconvolution DREAM Challenge](https://www.synapse.org/#!Synapse:syn15589870/wiki/), Deconer also provides a function to generate 'coarse' and 'fine' level mixture samples.

For the **'coarse'** level, Deconer generates mixture samples that consist of the following 8 cell types:

- B cells
- CD4 T cells
- CD8 T cells
- endothelial cells
- macrophages
- monocytes
- neutrophils
- NK cells

For the **'fine'** level, there will be 14 cell types:

- memory B cells
- naive B cells
- memory CD4 T cells
- active CD4 T cells
- regulatory T cells
- memory CD8 T cells
- naive CD8 T cells
- NK cells
- neutrophils
- monocytes
- myeloid dendritic cells
- macrophages
- fibroblasts
- endothelial cells

The following demonstration shows how to generate simulated samples. Since different deconvolution methods may be implemented in different programming languages, we provide the simulated data output as CSV files, which can be easily loaded into R, Python, as well as MATLAB.

```R
library(Deconer)

exprSim(n_sample = 50,  # generate 50 samples
        type = 'coarse',  # can be changed to 'fine'
        transform = 'TPM',
        outputPath = "./exprSim",
        mix_name = "coarse_mix.csv",  # mixture samples
        ref_name = "coarse_ref.csv",   # cell type specific reference
        prop_name = "coarse_prop.csv",  # groundtruth proportions
        refVar_name = "coarse_refVar.csv",  # expression variance for cell type specific reference
        train_name = "train_data.csv")  # data for generating reference, can be used for differential analysis
```

All the output files will be saved in the folder named 'exprSim'. The proportions will be randomly generated from a uniform distribution ([Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)). Additionally, Deconer provides output data for generating the external reference (parameter '**train_name**' in the function **exprSim**), which can be utilized for marker gene selection in differential expression analysis.

### Section 3.2: Gnerating Pseudo Bulk Data From scRNA-seq data

With the advancements in single-cell technologies, the prediction of cell type proportions has become more accurate. Numerous single-cell-based methods have been proposed, such as [MuSiC](https://doi.org/10.1038/s41467-018-08023-x) and [SCDC](https://doi.org/10.1093/bib/bbz166). Furthermore, the construction of in silico bulk data has become more straightforward.

For comprehensive simulations, Deconer provides two types of simulated bulk data.

The first type is the [human PBMC data from 10X](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k). This single-cell dataset contains more than 10,000 single cells from human blood. We processed this dataset using [10X Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) and [muon](https://muon.readthedocs.io/en/latest/). Cell type annotation was performed following the tutorial ['Processing gene expression of 10k PBMCs'](https://muon-tutorials.readthedocs.io/en/latest/single-cell-rna-atac/pbmc10k/1-Gene-Expression-Processing.html). We annotated 13 cell types as follows:

- intermediate_mono
- CD8+_naïve_T
- mDC
- CD4+_naïve_T
- NK
- memory_B
- CD4+_memory_T
- CD16_mono
- pDC
- naïve_B
- CD8+_activated_T
- CD14_mono
- MAIT

The second dataset is derived from mouse tissue. Tissue-level deconvolution is also an important aspect, such as predicting tissue origin from human liquid biopsy. [Han, *et al.*](https://doi.org/10.1016/j.cell.2018.02.001) constructed a high-quality scRNA-seq atlas for mice, which Deconer utilizes for in silico mixing. To ensure a robust simulation, Deconer incorporates data from 7 tissues of the female fetal mouse, with 1500 cells for each tissue.

- stomach
- lung
- liver
- kidney
- intestine
- brain
- gonad

The method for generating bulk data is similar with Section 3.1.

```R
scExprSim(
    n_sample = 50,  # generate 50 samples
    cell_number = 3000,
    ref_cell_number = 1000,
    p = 2/3,  # the proportion of cell number of building reference for each cell type
    transform = "TPM",
    outputPath = "./scExprSim",
    bulk_name = "scMouse_gene_expr.csv",
    ref_bulk_name = "scMouse_ref_bulk.csv",
    ref_sc_name = "scMouse_ref_sc.csv",
    ref_sc_label = "scMouse_ref_sc_label.csv",
    prop_name = "scMouse_prop.csv",
    type = "mouse_tissue"  ## 'mouse_tissue' or 'human_PBMC'
)
```

### Code Demo 1: Evaluating The Deconvolution Results In A Simple Manner

Generally, when a new deconvolution approach is proposed, it is essential to conduct comparisons among different methods. Deconer offers multi-level comparison and visualization functions. In this section, we will demonstrate how to perform a simple evaluation of deconvolution results.

For the demonstration, we have selected [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html), [DeconRNASeq](http://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html), and [FARDEEP](https://cran.r-project.org/web/packages/FARDEEP/index.html) from Bioconductor as the example methods. Users can follow the steps below to easily reproduce the results.

First, generate simulated expression data. In this example, we will use the massive RNA-seq dataset as a reference for generating the simulated data.

```R
library(Deconer)

# generate bulk data
exprSim(n_sample = 50,  # generate 50 samples
        type = 'coarse',  # can be changed to 'fine'
        transform = 'TPM',
        outputPath = "./exprSim",
        mix_name = "coarse_mix.csv",  # mixture samples
        ref_name = "coarse_ref.csv",   # cell type specific reference
        prop_name = "coarse_prop.csv",  # groundtruth proportions
        refVar_name = "coarse_refVar.csv",  # expression variance for cell type specific reference
        train_name = "train_data.csv")  # data for generating reference, can be used for differential analysis
```

Next, utilize [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html), [DeconRNASeq](http://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html), and [FARDEEP](https://cran.r-project.org/web/packages/FARDEEP/index.html) to perform deconvolution on the bulk data and save the results.

```R
library(EpiDISH)
library(DeconRNASeq)
library(FARDEEP)

mix <- read.csv(file = "./exprSim/coarse_mix.csv", header = T, row.names = 1)
ref <- read.csv(file = "./exprSim/coarse_ref.csv", header = T, row.names = 1)

# For a fast demo, we use only 500 markers.
mix <- as.matrix(mix[1:500, ])
ref <- as.matrix(ref[1:500, ])

# deconvolute with CIBERSORT algorithm
res1 <- epidish(beta.m = mix, ref.m = ref, method = "CBS")
p1 <- t(res1$estF)
# deconvolute with RPC algorithm
res2 <- epidish(beta.m = mix, ref.m = ref, method = "RPC")
p2 <- t(res2$estF)
# deconvolute with DeconRNASeq algorithm
res3 <- DeconRNASeq(datasets = as.data.frame(mix), signatures = as.data.frame(ref))
p3 <- t(res3$out.all)
colnames(p3) <- colnames(mix)
# deconvolute with FARDEEP algorithm
res4 <- fardeep(X = ref, Y = mix)
p4 <- t(res4$relative.beta)

# save the results
write.csv(x = p1, file = paste0("./exprSim/CBS.csv"), row.names = T, quote = F)
write.csv(x = p2, file = paste0("./exprSim/RPC.csv"), row.names = T, quote = F)
write.csv(x = p3, file = paste0("./exprSim/DeconRNASeq.csv"), row.names = T, quote = F)
write.csv(x = p4, file = paste0("./exprSim/FARDEEP.csv"), row.names = T, quote = F)
```

Next, users can compare the deconvolution results from various perspectives.

For instance, one can directly compare the root mean square error (RMSE) of each cell type.

```R
actual <- "./exprSim/coarse_prop.csv"
predicted <- c("./exprSim/CBS.csv",
               "./exprSim/RPC.csv",
               "./exprSim/DeconRNASeq.csv",
               "./exprSim/FARDEEP.csv")
label <- c("CBS", "RPC", "DeconRNASeq", "FARDEEP")

plot_multiple(actual = actual,
              predicted = predicted,
              label = label,
              method = "rmse",
              type = "celltype",  # can be changed to 'sample'
              figure = "boxplot")
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/box_ct_rmse.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Boxplot of rmse value for each cell type</div>
</center>

We can also compare the root mean square error (RMSE) of each sample by changing the parameter 'type' from 'celltype' to 'sample'.

Additionally, the overall RMSE can be visualized using a barplot.


```R
plot_multiple(actual = actual,
              predicted = predicted,
              label = label,
              method = "rmse",
              type = "all",
              figure = "barplot")
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/bat_all_rmse.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Barplot for overall rmse of each method</div>
</center>

Furthermore, the results can be visualized in a more specific manner using a heatmap.


```R
plot_multiple(actual = actual,
              predicted = predicted,
              label = label,
              method = "rmse",
              type = "celltype",
              figure = "heatmap")
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/heatmap_ct_rmse.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Heatmap for cell type rmse of each method</div>
</center>

Additionally, a scatter plot with the Pearson Correlation Coefficient (PCC) for each method can also be generated.

```R
plot_multiple(actual = actual,
              predicted = predicted,
              label = label,
              method = "pearson",
              type = "celltype",  # the color is assigned for each cell type
              figure = "scatterplot",
              nrow = 2)
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/scatter_ct_pearson.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Scatter plot for each method</div>
</center>

A combined heatmap with circles can be utilized to illustrate multiple metrics simultaneously.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/cheatmap_ct.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Heatmap with circles</div>
</center>

If users wish to plot the results for a single method, they can utilize the 'plot_single' function. For more detailed information, please refer to the [Deconer manual](https://honchkrow.github.io/Deconer/inst/documents/Deconer_manual.pdf).

## Section 4: Noise Analysis

In expression data analysis, both technical and biological noise cannot be overlooked. The presence of noise in bulk data can have a detrimental impact on deconvolution. To assess the stability and accuracy of different methods, Deconer offers functions to introduce noise with various structures into bulk data.

### Section 4.1: Generating Noised Bulk Data With Different Models

Various noise models have been used in different studies, including the normal, log-normal, and negative binomial models ([Jin, *et al.*](https://doi.org/10.1186/s13059-021-02290-6)). Deconer has implemented these models and offers a flexible interface for users to utilize.

- **Normal Model**

$$\textbf{M} = 2^{log_{2}(\textbf{R} \times \textbf{p} + 1) + N(0, \sigma p_t)}$$

- **Log-Nomal Model**

$$\textbf{M} = \textbf{R} \times \textbf{p} + 2^{N(0, \sigma p_t)}$$

In the aforementioned models, $\textbf{M}$ denotes the mixture sample, while $\textbf{R}$ and $\textbf{p}$ represent the reference and proportion, respectively. The level of noise is controlled by the product of a constant parameter $\sigma$ and a perturbation level parameter $p_t$ ([Jin, *et al.*](https://doi.org/10.1186/s13059-021-02290-6)). As indicated in the previous study, $\sigma$ is set to 10, and $p_t$ ranges from 0 to 1 in increments of 0.1.


- **Negative Binomial Model**

$$\mu_{i0} = r_{i0} \times L_{j}$$

$$\mu_{ij} = Gamma(shape = \frac{1}{\sigma_{i}^{2}}, scale = \frac{\mu_{i0}}{shape})$$

$$sigma_{i} = (1.8 \times p_{t} + \frac{1}{\sqrt{\mu_{i0}}}) \times exp(\frac{\delta}{2}) \qquad \delta \sim N(0, 0.25)$$

$$v_{ij} = Possion(\mu_{ij})$$

The simulation strategy was proposed by [Jin, *et al.*](https://doi.org/10.1186/s13059-021-02290-6) and [Law, *et al.*](https://doi.org/10.1186/gb-2014-15-2-r29). The noise level, denoted as $p_t$, can be controlled from 0 to 1. The expected genomic feature proportion of gene i in a cellular component is represented as $r_{i0}$. $L_{j}$ corresponds to the library size of sample $j$, and $\mu_{i0}$ represents the expected gene expression in the simulation. Two layers of variance are added using the gamma distribution and Poisson distribution. For more detailed information, please refer to [this article](https://doi.org/10.1186/s13059-021-02290-6).

The function "addNoiseExpr" can generate noise based on a normal, log-normal, or negative binomial model ([Jin, *et al.*](https://doi.org/10.1186/s13059-021-02290-6)). We highly recommend utilizing the **negative binomial** model. Below is a simple example code:

```R
# first, generate a simple dataset with 50 samples and 200 marker genes with function 'pseudoExpr'
res <- pseudoExpr(n_sample = 50, n_gene = 200)

# users can check the dimension of output
dim(res$mix)
dim(res$ref)
dim(res$prop)

# save the bulk data
write.csv(x = res$mix, file = "./addNoiseExpr/mix.csv", row.names = T, quote = F)
write.csv(x = res$ref, file = "./addNoiseExpr/ref.csv", row.names = T, quote = F)
write.csv(x = res$prop, file = "./addNoiseExpr/prop.csv", row.names = T, quote = F)

# the new data will be generated in a folder
addNoiseExpr(exprFile = "./addNoiseExpr/mix.csv",
             Pt = seq(0.1, 1, 0.1),  # parameter to control noise level
             type = "NB")  # "NB", "N" or "LN". 3 types of model.
```

A folder named "mix" will be generated in the path './addNoiseExpr', and bulk data with different levels of noise will be saved in separate files within the folder, as shown below.


```shell
mix.csv
ref.csv
prop.csv
mix/  
  ├── mix_NL_0.csv     # original data without in silico noise
  ├── mix_NL_0.1.csv   # noise level 0.1
  ├── mix_NL_0.2.csv 
  ├── mix_NL_0.3.csv 
  ├── mix_NL_0.4.csv 
  ├── mix_NL_0.5.csv 
  ├── mix_NL_0.6.csv 
  ├── mix_NL_0.7.csv 
  ├── mix_NL_0.8.csv 
  ├── mix_NL_0.9.csv 
  └── mix_NL_1.csv     # noise level 1
```

Afterward, users can easily test the performance of different methods. The following demonstration illustrates how to analyze the stability of different methods with varying noise levels.

### Code Demo 2: Evaluating The Deconvolution Results For Multiple Method

Comparing the deconvolution performance among different algorithms can help identify the most suitable method for a specific biological scenario.

In this example, we demonstrate how to utilize Deconer for conducting a systematic comparison of different methods using stability analysis. To provide a direct and efficient example, we have also employed [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html), [DeconRNAseq](https://www.bioconductor.org/packages/release/bioc/html/DeconRNASeq.html), and [FARDEEP](https://cran.r-project.org/web/packages/FARDEEP/index.html) for deconvolution.

First, perform deconvolution using different methods.

```R
library(EpiDISH)
library(DeconRNASeq)
library(FARDEEP)

ref <- as.matrix(read.csv(file = "./addNoiseExpr/ref.csv", header = T, row.names = 1))

for (NL in seq(0, 1, 0.1)) {
    writeLines(paste("Now, processing noise level:", NL, sep = " "))
    mixFile <- paste0("./addNoiseExpr/mix/mix_NL_", NL, ".csv")
    mix <- as.matrix(read.csv(file = mixFile, header = T, row.names = 1))

    # deconvolute with CIBERSORT algorithm
    res1 <- epidish(beta.m = mix, ref.m = ref, method = "CBS")
    p1 <- t(res1$estF)
    # deconvolute with RPC algorithm
    res2 <- epidish(beta.m = mix, ref.m = ref, method = "RPC")
    p2 <- t(res2$estF)
    # deconvolute with DeconRNASeq algorithm
    res3 <- DeconRNASeq(datasets = as.data.frame(mix), signatures = as.data.frame(ref))
    p3 <- t(res3$out.all)
    colnames(p3) <- colnames(p1)  # assign column name
    # deconvolute with FARDEEP algorithm
    res4 <- fardeep(X = ref, Y = mix)
    p4 <- t(res4$relative.beta)

    # save the results
    write.csv(x = p1, file = paste0("./addNoiseExpr/mix/CBS_", NL, ".csv"), row.names = T, quote = F)
    write.csv(x = p2, file = paste0("./addNoiseExpr/mix/RPC_", NL, ".csv"), row.names = T, quote = F)
    write.csv(x = p3, file = paste0("./addNoiseExpr/mix/DeconRNASeq_", NL, ".csv"), row.names = T, quote = F)
    write.csv(x = p4, file = paste0("./addNoiseExpr/mix/FARDEEP_", NL, ".csv"), row.names = T, quote = F)
}
```

First, we aim to assess the performance of a single method under varying levels of noise. Let's consider the CIBERSORT algorithm provided by [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html) as an example.

We can generate the trend of root mean square error (RMSE) across different noise levels. Additionally, we can plot cell type-specific metrics to visualize the <font color=red>estimation bias</font> for each <font color=red>cell type</font>.

```R
# the real proportion
actual <- "./addNoiseExpr/prop.csv"

# predicted proportions
predicted <- paste0("./addNoiseExpr/mix/CBS_", seq(0, 1, 0.1), ".csv")

# label for each file
noise_level <- paste0("NL", seq(0, 1, 0.1))

# boxplot
plot_multiple(actual = actual,
              predicted = predicted,
              label = noise_level,
              method = "rmse",
              type = "celltype",
              figure = "boxplot")

# heatmap
plot_multiple(actual = actual,
              predicted = predicted,
              label = noise_level,
              method = "rmse",
              type = "celltype",
              figure = "heatmap")
```

The output figures are presented below.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/noise_rmse_single_boxplot.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Box plot for rmse value of CIBERSORT</div>
</center>

It is evident that as the noise power increases, the quality of the deconvolution results deteriorates.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/noise_rmse_single_heatmap.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Cell type specific mape value for CIBERSORT</div>
</center>

Certainly, the results specific to cell types are influenced by noise power. **However, this issue has not been fully studied yet.**

Next, we will demonstrate how to perform a comprehensive cross-comparison among different methods.

```R
# the real proportion
actual <- "./addNoiseExpr/prop.csv"

# predicted proportions as a list
RPC <- paste0("./addNoiseExpr/mix/RPC_", seq(0, 1, 0.1), ".csv")
CBS <- paste0("./addNoiseExpr/mix/CBS_", seq(0, 1, 0.1), ".csv")
DeconRNASeq <- paste0("./addNoiseExpr/mix/DeconRNASeq_", seq(0, 1, 0.1), ".csv")
FARDEEP <- paste0("./addNoiseExpr/mix/FARDEEP_", seq(0, 1, 0.1), ".csv")
predicted <- list(RPC = RPC,
                  CBS = CBS,
                  DeconRNASeq = DeconRNASeq,
                  FARDEEP = FARDEEP)

noise_level <- paste0("NL_", seq(0, 1, 0.1))

# boxplot
plot_multiple2(actual = actual,
               predicted = predicted,
               condition = noise_level,
               method = 'rmse',  # compute rmse
               type = 'celltype',
               figure = 'boxplot')

# heatmap
plot_multiple2(actual = actual,
               predicted = predicted,
               condition = noise_level,
               method = 'mape',  # compute mape
               type = 'celltype',
               figure = 'heatmap')

```

The boxplot and heatmap for RMSE and MAPE are presented below.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/noise_rmse_multi_boxplot.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Boxplot of rmse for different deconvolution method</div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/noise_mape_multi_heatmap.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Heatmap of mape for different deconvolution method</div>
</center>

Typically, relying on a single metric may not be adequate to assess the efficacy of deconvolution, particularly for rare components ([Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)). Deconer offers a heatmap with circles that illustrates dual metrics in a single figure.

```R
# generate figure for rmse and pearson 
plot_multiple2(actual = actual,
               predicted = predicted,
               condition = noise_level,
               method = 'rmse',
               method2 = 'mape',
               type = 'celltype',
               figure = 'cheatmap')
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/noise_multi_cheatmap.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">circle heatmap of rmse and mape for different deconvolution method</div>
</center>

## Section 5: Rare Component Analysis

In the cell type deconvolution problem, cell types with extremely small proportions are often ignored. However, these cell types can play vital roles in certain situations, such as tumor-infiltrating lymphocytes (TILs), which exhibit low fractions in many cancer tissues. To address this issue, several methods have been proposed, such as [DWLS](https://doi.org/10.1038/s41467-019-10802-z) and [ARIC](https://doi.org/10.1093/bib/bbab362).

Deconer offers the functions 'rareExprSim' and 'rarescExprSim' for simulating bulk data with rare components. To conduct a comprehensive analysis, Deconer iterates over all cell types as potential rare components in a loop, using a pre-defined rare proportion gradient.

Here is an example involving 6 different deconvolution algorithms.

First, generate an in silico rare proportion dataset.

```R
# rare proportion is set to 0.001, 0.003, 0.005, 0.008, 0.01, 0.03 and 0.05
rareExprSim(p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05),
            outputPath = "./function_test/rareExprSim/",
            type = "coarse",
            transform = "TPM")
```

Next, we perform deconvolution on the simulated bulk data using different algorithms. Here, we present the deconvolution results obtained from [ARIC](https://doi.org/10.1093/bib/bbab362), [CIBERSORT](https://doi.org/10.1038/nmeth.3337), [EPIC](https://doi.org/10.1007/978-1-0716-0327-7_17), [dtangle](https://doi.org/10.1093/bioinformatics/bty926), [FARDEEP](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006976), and [DeconRNAseq](https://doi.org/10.1093/bioinformatics/btt090) as examples.

For simplicity, the specific code for each deconvolution method is not shown here. The following plotting code is based on known deconvolution results.

We can generate scatter plots for each method, as shown below.

```R
p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)

actual <- "coarse_prop.csv"

# use EPIC as an example
predicted <- "EPIC_coarse_coarse_p-0.01_fc-2.csv"

plot_rare(actual = actual,
          predicted = predicted,
          p_rare = p_rare,
          method = "spearman",
          celltype = TRUE,
          figure = "scatterplot")
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/scatter_R.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Scatter plot of rare component for EPIC</div>
</center>

From the above figure, users can analyze the deconvolution performance for each cell type.

Additionally, we can generate heatmaps and circular heatmaps to facilitate cross-comparisons.

```R
actual <- "./coarse_prop.csv"
predicted <- Sys.glob(paths = "./*_coarse_coarse_p-0.01_fc-2.csv")
label <- c("ARIC", "CBS", "EPIC", "dRNAseq", "dtangle", "fardeep")
p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)

plot_rare(actual = actual,
          predicted = predicted,
          p_rare = p_rare,
          method = "rmse",
          figure = "heatmap")

plot_rare(actual = actual,
          predicted = predicted,
          p_rare = p_rare,
          method = "rmse",
          method2 = "mape",
          figure = "cheatmap")
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/rare_multi_heatmap.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Heatmap plot for rmse</div>
</center>

<br />

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./inst/figures/rare_multi_cheatmap.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">circle heatmap plot for rmse and mape</div>
</center>

## Section 6: Single Cell Related Functions

Functions related to single-cell data are similar to those used for large-scale bulk data, such as 'scExprSim' and 'rarescExprSim'. For more detailed information, please consult the [Deconer manual](https://honchkrow.github.io/Deconer/inst/documents/Deconer_manual.pdf).

## Section 7: Well-Characterized Deconvolution Datasets

Additionally, we have curated 14 well-characterized deconvolution datasets for users. Some of these datasets include known cell-type proportions, while others may not have known true proportions but have associated phenotypic information. We provide these datasets in a processed format, including **bulk data**, **reference data**, and **true proportions**, making them readily usable. The summarized datasets are presented in the following table.


| NO. | NAME | Description | Reference Type | Proportion | Sample Size | Source |
| :----: | :----: | :----: | :----: | :----: | :----: | :----: |
| 1 | Abbas | Microarray | Bulk | known | 12 | [Abbas, *et al.*](https://doi.org/10.1371/journal.pone.0006098) |
| 2 | Becht | Microarray | Bulk | known | 10 | [Becht, *et al.*](https://doi.org/10.1186/s13059-016-1070-5) |
| 3 | Gong | Microarray | Bulk | known | 9 | [Gong, *et al.*](https://doi.org/10.1371/journal.pone.0027156) |
| 4 | Kuhn | Microarray | Bulk | known | 10 | [Kuhn, *et al.*](https://doi.org/10.1038/nmeth.1710) |
| 5 | Linsley | RNA-seq | Bulk | known | 5 | [Linsley, *et al.*](https://doi.org/10.1371/journal.pone.0109760) |
| 6 | Liu | RNA-seq | Bulk | known | 24 | [Liu, *et al.*](https://doi.org/10.1093/nar/gkv412) |
| 7 | Parsons | RNA-seq | Bulk | known | 30 | [Parsons, *et al.*](https://doi.org/10.1186/s12864-015-1912-7) |
| 8 | Shen-Orr | Microarray | Bulk | known | 33 | [Shen-Orr, *et al.*](https://doi.org/10.1038/nmeth.1439) |
| 9 | Shi | Microarray | Bulk | known | 60 | [Shi, *et al.*](https://doi.org/10.1038/nbt1239) |
| 10 | T2D | RNA-seq | Single Cell | unknown | 89 | [Fadista, *et al.*](https://doi.org/10.1073/pnas.1402665111) |
| 11 | TCGA_LUSC | RNA-seq | Bulk | unknown | 130 | [Vasaikar, *et al.*](https://doi.org/10.1093/nar/gkx1090) |
| 12 | TCGA_OV | RNA-seq | Bulk | unknown | 514 | [Vasaikar, *et al.*](https://doi.org/10.1093/nar/gkx1090) |
| 13 | kidney_Arvaniti | RNA-seq | Single Cell | unknown | 11 | [Arvaniti, *et al.*](https://doi.org/10.1038/srep26235) |
| 14 | kidney_Arvaniti_TPM | RNA-seq | Bulk | unknown | 11 | [Arvaniti, *et al.*](https://doi.org/10.1038/srep26235) |
| 15 | kidney_Craciun | RNA-seq | Single Cell | unknown | 19 | [Craciun, *et al.*](https://doi.org/10.1681/ASN.2015020225) |
| 16 | kidney_Craciun_TPM | RNA-seq | Bulk | unknown | 19 | [Craciun, *et al.*](https://doi.org/10.1681/ASN.2015020225) |
| 17 | TCGA 35 cancer datasets | RNA-seq | Bulk | unknown | - | [Vasaikar, *et al.*](https://doi.org/10.1093/nar/gkx1090) |

We have provided processed 'LUSC' and 'OV' datasets from TCGA, which have been validated by [FARDEEP](https://cran.r-project.org/web/packages/FARDEEP/index.html). Additionally, datasets related to other cancers are available in the form of expression matrix along with the original information.

Users can download these datasets from the following [page](https://github.com/Honchkrow/Deconer_dataset).

Note: Some datasets are collected from [dtangle](https://doi.org/10.1093/bioinformatics/bty926) and [MuSiC](https://doi.org/10.1038/s41467-018-08023-x).

Please cite the corresponding article when using these datasets.

## Citation

Zhang, Wei, et al. "Deconer: A comprehensive and systematic evaluation toolkit for reference-based cell type deconvolution algorithms using gene expression data." bioRxiv (2023): 2023-12.
