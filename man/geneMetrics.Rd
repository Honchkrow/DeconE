% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneMetrics.R
\name{geneMetrics}
\alias{geneMetrics}
\title{Compute SSIM for two expression matrix.}
\usage{
geneMetrics(expression_matrix, predicted_expression_matrix, method = NULL)
}
\arguments{
\item{expression_matrix}{Real expression data.
row: gene, column: samples.}

\item{predicted_expression_matrix}{The predicted proportion of cell types in matrix.
row: cell types, column: samples.}

\item{method}{One of c("SSIM", "JS"),}
}
\value{
vector
}
\description{
Compute SSIM for two expression matrix.
}
