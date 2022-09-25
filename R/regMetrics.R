#' @title Compute regression Metrics for matrix results.
#'
#' @description Compute several regression Metrics for cell type deconcolution.
#'
#' @param actual The groundcloth proportion of cell types in matrix.
#' row: cell types, column: samples.
#' @param predicted The predicted proportion of cell types in matrix.
#' row: cell types, column: samples.
#' @param method "rmse", "mape", "mae", "pearson", "spearman", "pearson2", "spearman2".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#' For pearson2 and spearman2, all data in matrix will be taken into
#' consideration, which mean only one number will be reported.
#'
#' @return a vector with the value for each sample.
#'
#' @importFrom Metrics rmse mape mae
#' @importFrom stats cor
#'
#' @export
#'
regMetrics <- function(actual, predicted, method = NULL) {
    if (!all(dim(actual) == dim(predicted))) {
        stop("Dimension of matrix actual and predicted are inconsistent.")
    } else {
        num_celltype <- nrow(actual)
        num_sample <- ncol(actual)
    }

    res <- c()
    if (method == "rmse") {
        for (i in seq(num_sample)) {
            res <- c(res, rmse(actual = actual[, i], predicted = predicted[, i]))
        }
        names(res) <- colnames(actual)
    } else if (method == "mape") {
        for (i in seq(num_sample)) {
            tmp_x <- actual[, i]
            tmp_y <- predicted[, i]
            idx <- which(tmp_x > 0)
            res <- c(res, mape(actual = tmp_x[idx], predicted = tmp_y[idx]))
        }
        names(res) <- colnames(actual)
    } else if (method == "mae"){
        for (i in seq(num_sample)) {
            res <- c(res, mae(actual = actual[, i], predicted = predicted[, i]))
        }
        names(res) <- colnames(actual)
    } else if (method == "pearson"){
        for (i in seq(num_sample)) {
            res <- c(res, cor(x = actual[, i], y = predicted[, i], method = "pearson"))
        }
        names(res) <- colnames(actual)
    } else if (method == "spearman") {
        for (i in seq(num_sample)) {
            res <- c(res, cor(x = actual[, i], y = predicted[, i], method = "spearman"))
        }
        names(res) <- colnames(actual)
    } else if (method == "pearson2"){
        res <- c(res, cor(x = as.vector(actual), y = as.vector(predicted), method = "pearson"))
    } else if (method == "spearman2") {
        res <- c(res, cor(x = as.vector(actual), y = as.vector(predicted), method = "spearman"))
    } else {
        stop("Parameter 'method' is invalid.")
    }
    return(res)
}
