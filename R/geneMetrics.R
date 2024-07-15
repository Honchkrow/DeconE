#' @title Compute SSIM for two expression matrix.
#'
#' @description Compute SSIM for two expression matrix.
#'
#' @param expression_matrix Real expression data.
#' row: gene, column: samples.
#' @param predicted_expression_matrix The predicted proportion of cell types in matrix.
#' row: cell types, column: samples.
#' @param method One of c("SSIM", "JS"),
#'
#' @return vector
#'
#' @export
#'
geneMetrics <- function(expression_matrix, predicted_expression_matrix, method = NULL) {
    # check matrix name
    if (check_matrix_names(expression_matrix, predicted_expression_matrix) == FALSE) {
        mess <- paste("Input error.")
        stop(mess)
    }

    # check dimension
    if (all(dim(expression_matrix) != dim(predicted_expression_matrix))) {
        mess <- paste("Input dimension error.")
        stop(mess)
    }

    n_genes <- nrow(expression_matrix)

    if (method == "SSIM") {
        scaled_expression_matrix <- expression_matrix / rowMaxs(expression_matrix)
        scaled_predicted_matrix <- predicted_expression_matrix / rowMaxs(predicted_expression_matrix)

        ssim_values <- numeric(n_genes)

        for (i in 1:n_genes) {
            ssim_values[i] <- compute_ssim(scaled_expression_matrix[i, ], scaled_predicted_matrix[i, ])
        }

        res <- ssim_values

    } else if (method == "JS") {
        prob_matrix <- calculate_probabilities(expression_matrix)
        predicted_prob_matrix <- calculate_probabilities(predicted_expression_matrix)

        js_values <- numeric(n_genes)
        for (i in 1:n_genes) {
            js_values[i] <- compute_js(prob_matrix[i, ], predicted_prob_matrix[i, ])
        }

        res <- js_values

    } else {
        mess <- paste("method error! must be one of the SSIM and JS.")
        stop(mess)
    }

    names(res) <- rownames(expression_matrix)

    return(res)

}


