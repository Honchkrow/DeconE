#' @title Plotting function for a single deconvolution method.
#'
#' @description Gnerating plots for a specific deconvolution method.
#'
#' @param actual The groundtruth proportion of cell types.
#' Can be a matrix or csv file with row and column names.
#' row: cell types, column: samples.
#' @param predicted The predicted proportion of cell types.
#' Can be a matrix or csv file with row and column names.
#' row: cell types, column: samples.
#' @param method One of the c("mae", "rmse", "mape", "kendall", "pearson", "spearman").
#' Note: for mape, cell types with real proportion 0 will be ignored.
#' Note: for scatter plot, method must in c("kendall", "pearson", "spearman").
#' @param type One of the c("sample", "celltpye", "all").
#' For "sample", generate metric for each sample. Scatter plot will assign different shape and color for each sample.
#' For "celltype", generate metric for each cell type. Scatter plot will assign different shape and color for each celltype.
#' For "all", generate metric for all data, which means flattening all data into an vector. Scatter plot will remove shape and color.
#' Note: boxplot is not supported for type="all".
#' @param figure One of the c("boxplot", "barplot", "scatterplot")
#' Note: boxplot is not supported for type="all".
#' Note: for scatter plot, method must in c("kendall", "pearson", "spearman").
#' @param errbar error bar type for barplot. One of the c("SD", "SE"), default: SE.
#' SD: Standard Deviation
#' SE: Standard Error
#'
#' @return The computed metrics as well as the plot data.
#'
#' @importFrom Metrics mape rmse mae
#' @importFrom ggplot2 ggplot geom_boxplot aes stat_boxplot theme
#' ggtitle xlab ylab element_text element_rect element_blank .data
#' geom_point scale_shape_manual geom_smooth labs geom_tile geom_text
#' coord_fixed scale_fill_gradient guides guide_legend guide_colorbar
#' guide_colourbar scale_x_discrete scale_radius theme_minimal position_dodge
#' geom_errorbar scale_fill_gradient2
#' @importFrom stats sd lm
#' @importFrom reshape2 melt
#' @importFrom ggpubr stat_cor
#'
#' @export
#'
#' @examples
#' res <- pseudoData(type = 1)
#' res <- plot_single(actual = res$actual, predicted = res$predicted, method = "mape", type = "sample", figure = "boxplot")
#'
plot_single <- function(actual,
                        predicted,
                        method = NULL,
                        type = "sample",
                        figure = "boxplot",
                        errbar = "SE") {
    p_true <- getInput(data = actual, name = "actual")
    p_pred <- getInput(data = predicted, name = "predicted")
    p_pred <- reorder_df(df1 = p_true, df2 = p_pred)


    if (method %in% c("mape", "mae", "rmse", 'smape', "kendall", "pearson", "spearman")) {
        res <- regMetrics(actual = p_true,
                          predicted = p_pred,
                          method = method,
                          type = type)
    } else {
        stop("Parameter method is invalid!")
    }

    if (figure == "boxplot") {
        if (type %in% c("sample", "celltype")) {
            title_mess <- paste("Boxplot for", type, method, sep = " ")
        } else {
            stop("Boxplot only supports type in c('sample', 'celltype')")
        }

        res <- data.frame(value = res)
        p <- ggplot(data = res, aes(x = "", y = .data$value)) +
            stat_boxplot(geom = "errorbar", width = 0.25) +
            geom_boxplot(fill = "#4271AE", outlier.colour = "red", outlier.shape = 19) +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
            ggtitle(title_mess) +
            xlab("") +
            ylab(method)
    } else if (figure == "barplot") {
        title_mess <- paste("Barplot for", type, method, sep = " ")
        # compute sd or se
        tmp_data <- sde(res, type = errbar)
        res <- data.frame(mean_value = mean(res), errbar_value = tmp_data)
        p <- ggplot(res) +
            geom_bar(aes(x = "", y = mean_value), stat = "identity", fill = "skyblue", alpha = 0.7) +
            geom_errorbar(aes(x = "", ymin = mean_value - errbar_value, ymax = mean_value + errbar_value), width = 0.4, colour = "orange", alpha = 0.9) +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
            ggtitle(title_mess) +
            xlab("") +
            ylab(method)


    } else if (figure == "scatterplot") {
        a1 <- melt(p_true)
        a2 <- melt(p_pred)
        data <- merge(x = a1, y = a2, by = c("Var1", "Var2"))
        colnames(data) <- c("celltype", "sample", "p_true", "p_pred")

        if (type == "sample") {
            p <- ggplot(data, aes(x = .data$p_true, y = .data$p_pred)) +
                geom_point(aes(shape = .data$sample, color = .data$sample)) +
                scale_shape_manual(values = seq(length(unique(data$sample)))) +
                geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
                stat_cor(method = method, p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
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
                     shape = "Samples", color = "Samples")
        } else if (type == "celltype") {
            p <- ggplot(data, aes(x = .data$p_true, y = .data$p_pred)) +
                geom_point(aes(shape = .data$celltype, color = .data$celltype)) +
                scale_shape_manual(values = seq(length(unique(data$celltype)))) +
                geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
                stat_cor(method = method, p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
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
            p <- ggplot(data, aes(x = .data$p_true, y = .data$p_pred)) +
                geom_point() +
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

    } else {
        stop("Parameter figure is invalid!")
    }

    return(list(data = res, plot = p))
}


















