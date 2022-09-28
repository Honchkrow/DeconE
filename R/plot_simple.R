#' @title Boxplot of deconvolution results for a specific method.
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
#' ggtitle xlab ylab element_text element_rect element_blank .data
#'
#' @export
#'
boxplot_simple <- function(actual,
                           predicted,
                           method = NULL) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = T,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (typeof(predicted) == "character") {
        p_pred <- as.matrix(read.csv(file = predicted,
                                     header = T,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (method %in% c("mape", "mae", "rmse", "pearson", "spearman")) {
        res <- regMetrics(actual = p_true, predicted = p_pred, method = method)
    } else {
        stop("Parameter must be mape, mae, rmse, pearson or spearman.")
    }

    res <- data.frame(value = res)

    p <- ggplot(data = res, aes(x = "", y = .data$value)) +
        # stat_boxplot(geom = "errorbar", width = 0.25) +
        geom_boxplot(fill = "#4271AE", outlier.colour = "red", outlier.shape = 19) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, size = 1),
              plot.title = element_text(size = 18, hjust = 0.5),
              axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
              axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5)) +
        ggtitle("Boxplot") +
        xlab("") +
        ylab(method)

    return(list(data = res, plot = p))
}





#' @title Scatter plot for proportion prediction.
#'
#' @description Generate scatter plot comparison for different cell types.
#'
#' @param actual The groundtruth proportion of cell types in matrix or csv file.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted The predicted proportion of cell types in matrix or csv file.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param method "pearson" (default), "kendall", or "spearman".
#' @param celltype TRUE or FALSE. Assign different type and color for different
#' cell types.
#'
#' @return The plot data.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_point scale_shape_manual geom_smooth
#' theme element_blank element_rect element_text labs .data
#' @importFrom stats lm
#' @importFrom ggpubr stat_cor
#'
#' @export
#'
scatter_simple <- function(actual, predicted, method, celltype = TRUE) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = T,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    } else {
        p_true <- actual
    }

    if (typeof(predicted) == "character") {
        p_pred <- as.matrix(read.csv(file = predicted,
                                     header = T,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    } else {
        p_pred <- predicted
    }

    if (!(method %in% c("pearson", "kendall", "spearman"))) {
        stop("Parameter must be pearson, kendall or spearman.")
    }

    a1 <- melt(p_true)
    a2 <- melt(p_pred)
    data <- merge(x = a1, y = a2, by = c("Var1", "Var2"))

    if (celltype) {
        p <- ggplot(data, aes(x = .data$value.x, y = .data$value.y)) +
            geom_point(aes(shape = .data$Var1, color = .data$Var1)) +
            scale_shape_manual(values = seq(length(unique(data$Var1)))) +
            geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
            stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  plot.title = element_text(size = 18, hjust = 0.5),
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 15)) +
            labs(title = "", x = "Ground Truth", y = "Prediction",
                 shape = "Cell Types", color = "Cell Types")
    } else {
        p <- ggplot(data, aes(x = .data$value.x, y = .data$value.y)) +
            geom_point() +
            scale_shape_manual(values = seq(length(unique(data$Var1)))) +
            geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
            stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  plot.title = element_text(size = 18, hjust = 0.5),
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 15)) +
            labs(title = "", x = "Ground Truth", y = "Prediction")
    }

    return(list(data = data, p = p))

}















