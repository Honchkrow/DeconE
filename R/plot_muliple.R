#' @title Plotting function for comparison multiple deconvolution method.
#'
#' @description Generating plots for multiple deconvolution method. This method
#' is designed for comparing the results from different methods under a certain
#' scenario or one method under different scenario. For example, comparing the
#' deconvolution effect of different methods from a specific noise level data, or
#' comparing the deconvolution effect of one method from different noise level.
#' Of course, there can be many samples in this certain scenario.
#' But for comparing the results from data with different noise
#' level, please use \code{\link{plot_multiple2}}.
#' Note: Function \code{\link{plot_multiple}} can reveal the deconvolution effect
#' for celltypes as well as samples directly. Function \code{\link{plot_multiple2}} reveals
#' the overall deconvolution results for different scenarios. The celltype or sample
#' specific effect can not be illustrated. Users should choose the appropriate
#' function. Of course, users can adopt the results from each function to
#' perform customized analysis.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a vector contains all the files with the predicted proportions
#' from different methods. Must be in csv file with row and column names.
#' row: cell types, column: samples.
#' For example, c("method1.csv", "method2.csv", "method3.csv").
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: base file name in parameter 'predicted'.
#' For example, c("method1", "method2", "method3").
#' @param method One of the c("mae", "rmse", "mape", "kendall", "pearson", "spearman").
#' Note: for mape, cell types with real proportion 0 will be ignored.
#' Note: for scatter plot, method must in c("kendall", "pearson", "spearman").
#' @param method2 One of the c("mae", "rmse", "mape", "kendall", "pearson", "spearman").
#' Generating the second metric in cheatmap.
#' Note: should be different from parameter method.
#' @param type One of the c("sample", "celltpye", "all").
#' For "sample", generate metric for each sample. Scatter plot will assign different shape and color for each sample.
#' For "celltype", generate metric for each cell type. Scatter plot will assign different shape and color for each celltype.
#' For "all", generate metric for all data, which means flattening all data into an vector. Scatter plot will remove shape and color.
#' Note: Parameter 'type' cannot be set to 'all' when plotting boxplot, heatmap and cheatmap.
#' @param figure One of the c("boxplot", "barplot", "scatterplot", "heatmap", "cheatmap")
#' Note: Parameter type cannot be set to 'all' when plotting boxplot, heatmap and cheatmap.
#' Note: for scatter plot, method must in c("kendall", "pearson", "spearman").
#' @param errbar error bar type for barplot. One of the c("SD", "SE"), default: SE.
#' SD: Standard Deviation
#' SE: Standard Error
#' @param nrow Only used in scatterplot. Control the layout of output figure.
#'
#' @return The computed metrics and plot data.
#'
#' @importFrom Metrics mape rmse mae
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.csv
#' @importFrom ggplot2 ggplot geom_boxplot aes stat_boxplot theme
#' ggtitle xlab ylab element_text element_rect element_blank .data
#' geom_point scale_shape_manual geom_smooth labs geom_tile geom_text
#' coord_fixed scale_fill_gradient guides guide_legend guide_colorbar
#' guide_colourbar scale_x_discrete scale_radius theme_minimal position_dodge
#' geom_errorbar
#' @importFrom stats sd lm
#' @importFrom reshape2 melt
#' @importFrom ggpubr stat_cor ggarrange
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 geom_bar scale_fill_gradient2
#'
#' @export
#'
plot_multiple <- function(actual,
                          predicted,
                          label = NULL,
                          method = NULL,
                          method2 = NULL,
                          type = "sample",
                          figure = "boxplot",
                          errbar = "SE",
                          nrow = 3) {

    check_method(method = method)

    if (is.null(label)) {
        label <- file_path_sans_ext(basename(predicted))
    }

    if (!(length(label) == length(predicted))) {
        stop("Different length between parameter 'predicted' and 'label'.")
    }

    # get real proportion
    p_true <- getInput(data = actual, name = "actual")

    # get predicted proportion
    res <- list()
    for (i in seq(length(predicted))) {
        this.label <- label[i]
        p_pred <- getInput(data = predicted[i], name = this.label)
        p_pred <- reorder_df(df1 = p_true, df2 = p_pred)

        res[[this.label]] <- regMetrics(actual = p_true,
                                        predicted = p_pred,
                                        method = method,
                                        type = type)
    }

    if (!(is.null(method2))) {
        res2 <- list()
        for (i in seq(length(predicted))) {
            this.label <- label[i]
            p_pred <- getInput(data = predicted[i], name = this.label)
            p_pred <- reorder_df(df1 = p_true, df2 = p_pred)

            res2[[this.label]] <- regMetrics(actual = p_true,
                                             predicted = p_pred,
                                             method = method2,
                                             type = type)
        }

    }

    if (figure == "boxplot") {
        if (type %in% c("sample", "celltype")) {
            title_mess <- paste("Boxplot for", type, method, sep = " ")
        } else {
            stop(" Parameter 'type' cannot be set to 'all' when plotting boxplot, heatmap and cheatmap. ")
        }
        df <- melt(data = do.call(cbind, res))
        colnames(df) <- c("type", "dmethod", "value")

        p <- ggplot(data = df, aes(x = .data$dmethod, y = .data$value)) +
            stat_boxplot(geom = "errorbar", width = 0.25) +
            geom_boxplot(width=0.25, outlier.shape = 19, fill = "#4271AE", outlier.colour = "red") +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  plot.title = element_text(size = 18, hjust = 0.5),
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5)) +
            ggtitle(title_mess) +
            xlab("Method") +
            ylab(method)

    } else if (figure == "barplot") {
        title_mess <- paste("Barplot for", type, method, sep = " ")
        if (type == "all") {
            df <- melt(res)
            colnames(df) <- c("value", "dmethod")
            p <- ggplot(df) +
                geom_bar(aes(x = dmethod, y = value), stat = "identity", fill = "skyblue", alpha = 0.7) +
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

        } else {
            data <- do.call(cbind, res)
            if (type == "sample") {
                mean_value <- apply(X = data, MARGIN = 1, FUN = mean)
                errbar_value <- apply(X = data, MARGIN = 1, FUN = Deconer:::sde, type = errbar)
            } else {
                mean_value <- apply(X = data, MARGIN = 2, FUN = mean)
                errbar_value <- apply(X = data, MARGIN = 2, FUN = Deconer:::sde, type = errbar)
            }

            df <- as.data.frame(cbind(mean_value, errbar_value))
            df <- rownames_to_column(df, "dmethod")

            p <- ggplot(df) +
                geom_bar(aes(x = dmethod, y = mean_value), stat = "identity", fill = "skyblue", alpha = 0.7) +
                geom_errorbar(aes(x = dmethod, ymin = mean_value - errbar_value, ymax = mean_value + errbar_value), width = 0.4, colour = "orange", alpha = 0.9) +
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

        }

    } else if (figure == "scatterplot") {
        # get predicted proportion
        p_list <- list()
        for (i in seq(length(predicted))) {
            this.label <- label[i]
            p_pred <- Deconer:::getInput(data = predicted[i], name = this.label)
            p_pred <- reorder_df(df1 = p_true, df2 = p_pred)

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

            p_list[[this.label]] <- p

        }

        df <- res
        p <- ggarrange(plotlist = p_list,
                       labels = label,
                       nrow = nrow,
                       ncol = ceiling(length(p_list)/nrow))



    } else if (figure == "heatmap") {
        if (!(type %in% c("sample", "celltype"))) {
            stop(" Parameter 'type' cannot be set to 'all' when plotting boxplot, heatmap and cheatmap. ")
        }

        df <- do.call(rbind, res)
        df <- round(x = df, digits = 4)
        df <- melt(df)

        if (type == "sample") {
            xx <- "Sample"
        } else {
            xx <- "Cell Type"
        }

        p <- ggplot(df, aes(x = .data$Var2, y = .data$Var1, fill = value)) +
            geom_tile(color = "black") +
            geom_text(aes(label = value), color = "white") +
            coord_fixed() +
            scale_fill_gradient(low = "blue", high = "red") +
            guides(fill = guide_colourbar(title = method)) +
            theme(plot.title = element_text(size = 18, hjust = 0.5),
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 15)) +
            labs(title = "", x = type, y = "Method")

    } else if (figure == "cheatmap") {
        if (!(type %in% c("sample", "celltype"))) {
            stop(" Parameter 'type' cannot be set to 'all' when plotting boxplot, heatmap and cheatmap. ")
        }

        fi1 <- do.call(rbind, res)
        fi2 <- do.call(rbind, res2)
        a1 <- melt(fi1)
        a2 <- melt(fi2)

        data <- merge(x = a1, y = a2, by = c("Var1", "Var2"))

        if (method %in% c("pearson", "spearman")) {
            trend_label <- c("Bad", "Good")
        } else {
            trend_label <- c("Good", "Bad")
        }

        df <- data
        p <- ggplot(data, aes(x = .data$Var2,
                              y = forcats::fct_rev(.data$Var1),
                              fill = .data$value.x,
                              size = .data$value.y)) +
            geom_point(shape = 21,
                       stroke = 0) +
            scale_x_discrete(position = "bottom") +
            scale_radius(range = c(5, 15)) +
            # scale_fill_gradient(low = "yellow",
            #                     high = "purple3",
            #                     breaks = c(min(data$value.x), max(data$value.x)),
            #                     # labels = trend_label,
            #                     labels = c(round(x = min(data$value.x), digits = 2),
            #                                round(x = max(data$value.x), digits = 2)),
            #                     limits = c(min(data$value.x), max(data$value.x))) +
            scale_fill_gradient2(low = "yellow",
                                 mid = "red",
                                 high = "mediumblue",
                                 midpoint = mean(data$value.x),
                                 breaks = c(min(data$value.x), max(data$value.x)),
                                 labels = c(round(x = min(data$value.x), digits = 2),
                                            round(x = max(data$value.x), digits = 2)),
                                 limits = c(min(data$value.x), max(data$value.x))) +
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
            guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25),
                                       label.position = "bottom",
                                       title.position = "right",
                                       order = 1),
                   fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
            labs(size = method2,
                 fill = method,
                 x = NULL,
                 y = NULL)

    } else {
        stop("Parameter figure is invalid!")
    }

    return(list(data = df, plot = p))

}
