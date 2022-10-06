# **decone**: An easy-to-use and comprehensive evaluation toolkit for cell type deconvolution from expression data

## Section 1: Introduction
Quantifying cell proportions is an important problem in bioinformatics. cell type proportion is related with certain phenotypes or diseases ([Wei, *et al.*](https://doi.org/10.1093/bib/bbab362)).

Here, we proposed a cell type <u>decon</u>volution <u>e</u>valuating toolkit named '**decone**' to perform comprehensive and systematic analysis for different algorithms.

**decone** consists of 5 main part functions as below.
- Pseudo bulk data generation (including bulk and single cell).
- Stability analysis under different types of noise.
- Rare component analysis.
- Unknown component analysis.
- Comprehensive evaluation metrics.

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


