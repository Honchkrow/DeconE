#' @title Plot circle heatmap for rare components.
#'
#' @description Boxplot for deconvolution results with multiple samples.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted The predicted proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param method "rmse", "mape", "mae", "pearson" or "spearman".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#'
#' @return The computed metrics and plot data.
#'
#' @importFrom Metrics mape rmse mae
#' @importFrom ggplot2 ggplot geom_boxplot aes stat_boxplot theme
#' ggtitle xlab ylab element_text element_rect element_blank
#'
#' @export
#'
