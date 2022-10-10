# **decone**: An easy-to-use and comprehensive evaluation toolkit for cell type deconvolution from expression data


**Links**:

- [decone documentaion](123)
- [decone vignettes](123)


## Section 1: Introduction
Cell type proportion is related with certain phenotype or disease ([Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)). Therefore, quantifying cell or tissue proportions is an important problem in bioinformatics.

Here, we proposed a cell type <u>decon</u>volution <u>e</u>valuating toolkit named '**decone**' to perform comprehensive and systematic analysis for different algorithms.

**decone** consists of 6 main part functions as below.
- Pseudo bulk data generation (including bulk and single cell).
- Stability analysis under different types of noise.
- Rare component analysis.
- Unknown component analysis.
- Comprehensive evaluation metrics.
- Well characterized datasets for deconvolution utilities.

In the following parts, we will introduce each function along with how to compute the evaluation metrics for comparison different deconvolution methods.

## Section 2: Installation

decone is based on R and can be eaisly installed in windows, linux as well as mac OS.

First, users should install [R](https://www.r-project.org/).

Next, install devtools and decone.

```
# install devtools
install.packages('devtools')

# install the decone package
devtools::install_github('Honchkrow/decone')

# load decone
library(decone)
```

## Section 3: Gnerating Pseudo Bulk Data

Generating pseudo bulk data is a challenging problem. Inspired by the former work ([Francisco, *et al.*](https://doi.org/10.1038/s41467-020-20288-9), [Wang, *et al.*](https://doi.org/10.1038/s41467-018-08023-x), [Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)) and [Tumor Deconvolution DREAM Challenge](https://www.synapse.org/#!Synapse:syn15589870/wiki/), decone provides different pseudo data generation strategies from bulk RNA-seq data as well as scRNA-seq data.

### Section 3.1: Gnerating Pseudo Bulk Data From Massive RNA-seq Studies

Many deconvolution methods need cell type specific bulk data as the prior knowledge during deconvolution. In order to have a more realistc simulation, we collected 302 well-characterized bulk RNA-seq data to generate pseudo bulk data. When generating data, 1/3 samples will be used for generating the mixture and the rest will be used for generating external reference.

Inspired by [Tumor Deconvolution DREAM Challenge](https://www.synapse.org/#!Synapse:syn15589870/wiki/), decone also provides to generate 'coarse' and 'fine' level mixture sample. 

For **'coarse'** level sample, decone generates mixture sample which contains the following 8 cell types.
- B cells
- CD4 T cells
- CD8 T cells
- endothelial cells
- macrophages
- monocytes
- neutrophils
- NK cells. 

For **'fine'** level, there will be 14 cell types.
- memory B cells
- naive B cells
- memory CD4 T cells
- aive CD4 T cells
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

The following demo shows how to generate simulated sample. Considering that methods may be written in different language, we output the simulated data as csv file.

```R
# create the folder 'test' before running
exprSim(n_sample = 50,  # generating 50 samples
        type = 'coarse',  # can be changed to 'fine'
        transform = 'TPM', 
        outputPath = "./test",
        mix_name = "coarse_mix.csv",  # simulated 50 mixture samples
        ref_name = "coarse_ref.csv",   # cell type specific reference
        prop_name = "coarse_prop.csv",  # simulated porportions
        refVar_name = "coarse_refVar.csv",  # expression variance for cell type specific reference
        train_name = "train_data.csv")  # data for generating reference

```

The proportion will be generated randomly from uniform distribution ([Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)). decone also output the data for generating the external reference which can be used for differential expression analysis in marker gene selection.


### Section 3.2: Gnerating Pseudo Bulk Data From scRNA-seq data

With the development of single cell technologies, cell type deconvolution could be more accurate. Lots of single cell based methods had been proposed such as [MuSiC](https://doi.org/10.1038/s41467-018-08023-x) and [SCDC](https://doi.org/10.1093/bib/bbz166). In addition, constructing in silico bulk data also becomes directly. [Han, *et al.*](https://doi.org/10.1016/j.cell.2018.02.001) built an scRNA-seq atlas for mouse with high quality and decone adopted this study to perform in silico mixing. For an solid simulation, decone adpoted 7 tissues from femal fetal mouse and 1500 cells for each, including stomach, lung, liver, kidney, intestine, brain and gonad.

The method for generating bulk data is similar with Section 3.1.

```R
# create the folder 'test' before running
scExprSim(n_sample = 50,
          p = 2 / 3,
          transform = "TPM",
          outputPath = "./test",
          mix_name = "scMouse_gene_expr.csv",
          ref_name = "scMouse_ref.csv",
          prop_name = "scMouse_prop.csv",
          train_name = "scMouse_ref_rawCount.csv")  # data for generating reference

```

### Section 3.3: Evaluating The Deconvolution Results For Single Method

In this part, we will demonstrate how to analysis the results for a single deconvolution method.

We adpoted [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html) which provides 3 types of deconvolution algorithms as demo.

First, generating simulated expression data.

```R
# install EpiDISH
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EpiDISH")

# generate bulk data
exprSim(n_sample = 50,  # generating 50 samples
        type = 'coarse',  # can be changed to 'fine'
        transform = 'TPM', 
        outputPath = "./test",
        mix_name = "coarse_mix.csv",  # simulated 50 mixture samples
        ref_name = "coarse_ref.csv",   # cell type specific reference
        prop_name = "coarse_prop.csv",  # simulated porportions
        refVar_name = "coarse_refVar.csv",  # expression variance for cell type specific reference
        train_name = "train_data.csv")  # data for generating reference
```

Second, use [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html) to deconvolute the bulk data and save the result.

```R
library(EpiDISH)

mix <- read.csv(file = "./test1/coarse_mix.csv", header = T, row.names = 1)
ref <- read.csv(file = "./test1/coarse_ref.csv", header = T, row.names = 1)

# For a fast demo, we use only 500 markers.
mix <- as.matrix(mix[1:500, ])
ref <- as.matrix(ref[1:500, ])

res <- epidish(beta.m = mix, ref.m = ref, method = "CBS")
prop_pred <- t(res$estF)

write.csv(x = prop_pred, file = "./test1/prop_pred.csv", row.names = T, quote = F)
```

Next, using decone to Evaluate the deconvolution results.

```R
actual <- "./test1/coarse_prop.csv"
predicted <- "./test1/prop_pred.csv"

# Plot boxplot for different metrics
# Evaluation method can be changed to "rmse", "mape", "mae", "pearson" or "spearman".
res <- boxplot_simple(actual = actual,
                      predicted = predicted,
                      method = "rmse")
res$plot
```

The rmse boxplot of 50 samples will be like below.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./inst/figures/boxplotSimple.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Boxplot of rmse value for 50 samples</div>
</center>

For generating the comparison plot between different method under a certain dataset, please check the **Section 4**.

Also, users can generate pretty scatter plot for each dataset like below.

```R
# Plot scatter-plot for different metrics
# Evaluation method can be changed to "pearson" (default), "kendall", or "spearman".
res <- scatter_simple(actual = actual,
                      predicted = predicted,
                      method = "spearman")
res$plot
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./inst/figures/scatterSimple.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Scatter-plot of rmse value for 50 samples</div>
</center>

Generally, when a new deconvolution approch is proposed, cross-comparison between different methods is needed. decone provides multi-level comparison and visualization functions. The cross-comparison function will be introduced along with the next parts.

## Section 4: Noise Analysis

In expression data analysis, technical and biological noise cannot be ignored. Noised existed in bluk data brings negative influence for deconvolution. In order to measure the stability and accuracy for different methods, decone provides functions to add noise with different structure on bulk data.


### Section 4.1: Generating Noised Bulk Data With Negative Binomial Model

Many the studies adopted normal or log-normal model to generate noised bulk data. However, as [Jin, *et al.*](https://doi.org/10.1186/s13059-021-02290-6) pointed, the negative binomial model recapitulates noise structures of real bulk data. 

Function "addNoiseExpr" can produce noise based on normal, log-normal or negative binomial model ([Jin, *et al.*](https://doi.org/10.1186/s13059-021-02290-6)). We strongly recommended using 'negative binomial' model. An simple example is provided as the following code.

```R
# generate a simulated dataset with 50 samples and 200 marker genes with function 'pseudoExpr'
res <- pseudoExpr(n_sample = 50, n_gene = 200)

# users can check the dimension of output
dim(res$mix)
dim(res$ref)
dim(res$prop)

# save the bulk data
write.csv(x = res$mix, file = "mix.csv", row.names = T, quote = F)
write.csv(x = res$ref, file = "ref.csv", row.names = T, quote = F)
write.csv(x = res$prop, file = "prop.csv", row.names = T, quote = F)

# the new data will be generated in a folder
addNoiseExpr(exprFile = "mix.csv", 
             Pt = seq(0.1, 1, 0.1),  # parameter to control noise level
             type = "NB")  # "NB", "N" or "LN". 3 types of model.
```

A folder named "mix" will be generated and bulk data with different level noise will be save in each file respectively like below.

``` 
mix.csv
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

After this, users can test the performance of different methods easily.


### Section 4.2: Evaluating The Deconvolution Results For Multiple Method

Comparison the deconvolution performance between different algorithms helps to find the appropriate method for a certain biology scenario.

Here, taken stability analysis as an example, we shows how to use decone to perform cross-comparison between different methods. In order to give a direct and fast example, we use [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html), [DeconRNAseq](https://www.bioconductor.org/packages/release/bioc/html/DeconRNASeq.html) as well as [FARDEEP](https://cran.r-project.org/web/packages/FARDEEP/index.html) to perform deconvolution.

```R
library(EpiDISH)
library(DeconRNASeq)
library(FARDEEP)

ref <- as.matrix(read.csv(file = "ref.csv", header = T, row.names = 1))

for (NL in seq(0, 1, 0.1)) {
    writeLines(paste("Now, processing noise level:", NL, sep = " "))
    mixFile <- paste0("./mix/mix_NL_", NL, ".csv")
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
    # deconvolute with FARDEEP algorithm
    res4 <- fardeep(X = ref, Y = mix)
    p4 <- t(res4$relative.beta)

    # save the results
    write.csv(x = p1, file = paste0("./mix/CBS_", NL, ".csv"), row.names = T, quote = F)
    write.csv(x = p2, file = paste0("./mix/RPC_", NL, ".csv"), row.names = T, quote = F)
    write.csv(x = p3, file = paste0("./mix/DeconRNASeq_", NL, ".csv"), row.names = T, quote = F)
    write.csv(x = p4, file = paste0("./mix/FARDEEP_", NL, ".csv"), row.names = T, quote = F)
}
```

First, we want to see the performance of a single method under different level noise. Let's take CIBERSORT algorithm provided by [EpiDISH](https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html) as an example.

We can generate the rmse trend along with different noise level. Also, we can plot the cell type specific metrics to see the <font color=red>estimation bias</font> for each cell types.

```R
# the real proportion
actual <- "prop.csv"

# predicted proportions 
predicted <- paste0("./mix/CBS_", seq(0, 1, 0.1), ".csv")

# generate box plot of rmse value
res <- boxplot_NGrad(actual = actual,
                     predicted = predicted,
                     label = paste0("NL_", seq(0, 1, 0.1)),
                     method = "rmse")

res$plot

# generate heatmap of mape for each cell types
res <- heatmap_NGradCT(actual = actual,
                       predicted = predicted,
                       label = paste0("NL_", seq(0, 1, 0.1)),
                       method = "mape")
res$plot
```

The output figures are as follows.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./inst/figures/boxplot_NGrad.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Box plot for rmse value of CIBERSORT</div>
</center>


It is clear that with the growth of noise power, the deconvolution results becomes worse.


<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./inst/figures/heatmap_NGradCT.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Cell type specific mape value for CIBERSORT</div>
</center>

Of course, cell types specific results are influenced by noise power. **However, this issue is still not fully studied.**

Next, we wan to have a comprehensive cross-comparison between different methods.

```R
# the real proportion
actual <- "prop.csv"

# predicted proportions as a list
RPC <- paste0("./mix/RPC_", seq(0, 1, 0.1), ".csv")
CBS <- paste0("./mix/CBS_", seq(0, 1, 0.1), ".csv")
DeconRNASeq <- paste0("./mix/DeconRNASeq_", seq(0, 1, 0.1), ".csv")
FARDEEP <- paste0("./mix/FARDEEP_", seq(0, 1, 0.1), ".csv")
predicted <- list(RPC = RPC,
                  CBS = CBS,
                  DeconRNASeq = DeconRNASeq,
                  FARDEEP = FARDEEP)

noise_level <- paste0("NL_", seq(0, 1, 0.1))

# boxplot
boxplot_NcrossCompare(actual = actual,
                      predicted = predicted,
                      label = noise_level,
                      method = "rmse")

# heatmap
heatmap_NcrossCompare(actual = actual,
                      predicted = predicted,
                      label = noise_level,
                      method = "mape")
```

The bosplot and heatmap for rmse and mape are as follows.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./inst/figures/boxplot_NcrossCompare.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Boxplot of rmse for different deconvolution method</div>
</center>


<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./inst/figures/heatmap_NcrossCompare.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Heatmap of mape for different deconvolution method</div>
</center>


Usually, only one metrics may not be sufficient to reveal the deconvolution efficacy ([Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)). decone provides circle heatmap for illustrating the double metrics in one figure.


```R
# generate figure for rmse and pearson 
cheatmap_NcrossCompare(actual = actual,
                       predicted = predicted,
                       label = noise_level,
                       method1 = "rmse",
                       method2 = "pearson")
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./inst/figures/cheatmap_NcrossCompare.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">circle heatmap of rmse and PCC for different deconvolution method</div>
</center>


## Section 5: Rare Component Analysis





## Section 5: Unknown Component Analysis



## Section 5: Well-Characterized Deconvolution Datasets




## References


























