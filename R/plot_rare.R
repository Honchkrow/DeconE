#' @title Scatter plot for rare component.
#'
#' @description Generate scatter plot for rare component.
#' Note: only the rare proportion will be considered.
#' The 'scatterplot' function can be used to estimate the deconvolution power for each cell type.
#' The 'heatmap' and 'cheatmap' functions can be used to compare the deconvolution power for each method.
#'
#' @param actual The groundtruth proportion of cell types in matrix or csv file.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted The predicted proportion of cell types in matrix or csv file.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' Note: For parameter 'figure' is scatterplot, 'predicted' must be one matrix or filename.
#' This means scatterplot is designed for generating plot for one method.
#' For example, predicted = 'method1_predicted.csv'
#' Note: For parameter 'figure' is in c("heatmap", "cheatmap"), 'predicted' must be a vector of filename.
#' This means heatmap is designed for generating plot for multiple methods.
#' For example, predicted = c('method1_predicted.csv', 'method2_predicted.csv', 'method3_predicted.csv')
#' @param label a vector contains the label corresponding to the predicted
#' proportion file from different methods. Default: base file name in parameter 'predicted'.
#' For example, c("method1", "method2", "method3").
#' Note: this parameter is only supported by the heatmap and cheatmap functions.
#' @param p_rare A vector of proportions. should be the same with function \code{\link{rareExprSim}}.
#' Default: c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)
#' @param method One of the c("mae", "rmse", "mape", "kendall", "pearson", "spearman").
#' For parameter 'figure' is scatterplot, method must be one of c("kendall", "pearson", "spearman")
#' Note: kendall, pearson and spearman is not supported by heatmap and cheatmap.
#' @param method2 One of the c("mae", "rmse", "mape", "kendall", "pearson", "spearman").
#' Generating the second metric in cheatmap.
#' Note: method2 will be disabled when figure is in c("scatterplot", "heatmap")
#' Note: kendall, pearson and spearman is not supported by heatmap and cheatmap.
#' @param type Not used. All metrics will be computed for each rare component,
#' regardless of sample and cell type.
#' @param celltype TRUE or FALSE. Assign different type and color for different
#' cell types. Default: TRUE.
#' Note: This parameter only supported by the scatterplot.
#' @param figure One of the c("scatterplot", "heatmap", "cheatmap").
#' For different type of figure, users should pass the correct parameters.
#'
#' @return The plot data.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_boxplot aes stat_boxplot theme
#' ggtitle xlab ylab element_text element_rect element_blank .data
#' geom_point scale_shape_manual geom_smooth labs geom_tile geom_text
#' coord_fixed scale_fill_gradient guides guide_legend guide_colorbar
#' guide_colourbar scale_x_discrete scale_radius theme_minimal position_dodge
#' geom_errorbar
#' @importFrom ggpubr stat_cor
#' @importFrom stats lm
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
#' @examples
#' rareExprSim()
#' # after the simulation, just pass the predicted file
#' # to the parameter "actual" and "predicted"
#'
plot_rare <- function(actual,
                      predicted,
                      p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05),
                      method = NULL,
                      method2 = NULL,
                      type = NULL,
                      celltype = TRUE,
                      figure = NULL) {
    if (!(method %in% c("mape", "mae", "rmse", "pearson", "spearman", "kendall"))) {
        stop("Parameter method must be one of c('mae', 'mape', 'rmse', 'pearson', 'spearman', 'kendall').")
    }

    # get real proportion
    p_true <- decone:::getInput(data = actual, name = "actual")

    if (figure == "scatterplot") {
        p_pred <- decone:::getInput(data = predicted, name = "predicted")

        p_len <- length(p_rare)

        pt_list <- list()
        pp_list <- list()

        for (i in seq(nrow(p_true))) {
            ct <- rownames(p_true)[i]
            idx_start <- (i - 1) * p_len + 1
            idx_end <- i * p_len
            pt_list[[ct]] <- p_true[ct, idx_start:idx_end]
            pp_list[[ct]] <- p_pred[ct, idx_start:idx_end]
        }

        pt_m <- do.call(rbind, pt_list)
        colnames(pt_m) <- p_rare
        pp_m <- do.call(rbind, pp_list)
        colnames(pp_m) <- p_rare

        df <- melt(pp_m)

        if (celltype) {
            p <- ggplot(df, aes(x = .data$Var2, y = .data$value)) +
                geom_point(aes(shape = .data$Var1, color = .data$Var1)) +
                scale_shape_manual(values = seq(length(unique(df$Var1)))) +
                geom_smooth(method = lm , color = "red", fill = "#69b3a2", se = TRUE) +
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
            p <- ggplot(df, aes(x = .data$Var2, y = .data$value)) +
                geom_point() +
                scale_shape_manual(values = seq(length(unique(df$Var1)))) +
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
                labs(title = "", x = "Ground Truth", y = "Prediction")
        }
    } else if (figure %in% c("heatmap", "cheatmap")) {
        if ((method %in% c("pearson", "spearman", "kendall"))) {
            stop(" kendall, pearson and spearman is not supported by heatmap and cheatmap. ")
        }


        if (is.null(label)) {
            label <- file_path_sans_ext(basename(predicted))
        }

        if (!(length(label) == length(predicted))) {
            stop("Different length between parameter 'predicted' and 'label'.")
        }

        p_len <- length(p_rare)

        res <- list()
        res2 <- list()
        for (i in seq(length(predicted))) {
            p_pred <- decone:::getInput(data = predicted[i], name = "predicted")
            tmp_label <- label[i]

            if (!all(dim(p_pred) == dim(p_true))) {
                mess <- paste("The dimension of input file",
                              predicted[1],
                              "is not consistent.")
                stop(mess)
            }

            pt_list <- list()
            pp_list <- list()

            # construct the matrix, row: celltype, column: rare proportion
            for (i in seq(nrow(p_true))) {
                ct <- rownames(p_true)[i]
                idx_start <- (i - 1) * p_len + 1
                idx_end <- i * p_len
                pt_list[[ct]] <- p_true[ct, idx_start:idx_end]
                pp_list[[ct]] <- p_pred[ct, idx_start:idx_end]
            }

            pt_m <- do.call(rbind, pt_list)
            colnames(pt_m) <- p_rare
            pp_m <- do.call(rbind, pp_list)
            colnames(pp_m) <- p_rare

            res[[tmp_label]] <- regMetrics(actual = pt_m,
                                           predicted = pp_m,
                                           method = method, type = "sample")
            if (!(is.null(method2))) {
                if ((method2 %in% c("pearson", "spearman", "kendall"))) {
                    stop(" kendall, pearson and spearman is not supported by heatmap and cheatmap. ")
                }
                res2[[tmp_label]] <- regMetrics(actual = pt_m,
                                                predicted = pp_m,
                                                method = method2, type = "sample")
            }

        }

        if (figure == "heatmap") {
            df <- as.matrix(do.call(cbind, res))
            df <- round(x = df, digits = 4)
            df <- melt(df)
            colnames(df) <- c('rare', 'Method', 'value')
            df$rare <- factor(df$rare)

            title_mess <- method

            p <- ggplot(df, aes(x = .data$Method, y = .data$rare, fill = value)) +
                geom_tile(color = "black") +
                geom_text(aes(label = value), color = "white") +
                coord_fixed() +
                scale_fill_gradient(low = "blue", high = "red") +
                guides(fill = guide_colourbar(title = title_mess))
        } else if (figure == "cheatmap") {
            df1 <- as.matrix(do.call(cbind, res))
            df1 <- round(x = df1, digits = 4)
            df1 <- melt(df1)
            colnames(df1) <- c('rare', 'Method', 'value1')
            df1$rare <- factor(df1$rare)

            df2 <- as.matrix(do.call(cbind, res2))
            df2 <- round(x = df2, digits = 4)
            df2 <- melt(df2)
            colnames(df2) <- c('rare', 'Method', 'value2')
            df2$rare <- factor(df2$rare)

            df <- merge(x = df1, y = df2, by = c("rare", "Method"))

            size_lab <- method2
            fill_lab <- method
            p <- ggplot(df, aes(x = .data$Method,
                                y = forcats::fct_rev(.data$rare),
                                fill = .data$value1,
                                size = .data$value2)) +
                geom_point(shape = 21,
                           stroke = 0) +
                scale_x_discrete(position = "bottom") +
                scale_radius(range = c(5, 15)) +
                scale_fill_gradient(low = "orange",
                                    high = "blue",
                                    breaks = c(min(df$value1), max(df$value1)),
                                    # labels = trend_label,
                                    labels = c(round(x = min(df$value1), digits = 2),
                                               round(x = max(df$value1), digits = 2)),
                                    limits = c(min(df$value1), max(df$value1))) +
                theme_minimal() +
                theme(panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.position = "bottom",
                      legend.text = element_text(size = 8),
                      legend.title = element_text(size = 8),
                      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
                      axis.text.y = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
                      panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                guides(size = guide_legend(override.aes = list(fill = "black", color = "black", stroke = 0.25),
                                           label.position = "bottom",
                                           title.position = "right",
                                           order = 1),
                       fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
                labs(size = size_lab,
                     fill = fill_lab,
                     x = NULL,
                     y = NULL)

        } else {
            stop(" Parameter 'figure' is not valid. ")
        }
    } else {
        stop(" Parameter 'figure' is not valid. ")
    }



    return(list(data = df, plot = p))

}





