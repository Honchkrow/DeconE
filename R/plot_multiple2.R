#' @title Plotting function for comparison multiple deconvolution method.
#'
#' @description Generating plots for multiple deconvolution method. This method
#' is designed for comparing the results from different methods under different
#' scenario. For example, comparing the deconvolution effect of different methods
#' from various noise level data. Of course, there can be many samples in a
#' certain scenario. But for comparing the results from data with a certain noise
#' level, please use \code{\link{plot_multiple}}.
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
#' @param predicted a list of vectors. Each vector contains all the files with
#' the predicted proportions in different scenarios. Must be in csv file with
#' row and column names.
#' row: cell types, column: samples.
#' For example, list(method1 = c("m1_noise0.csv", "m1_noise1.csv", "m1_noise2.csv"),
#' method2 = c("m2_noise0.csv", "m2_noise1.csv", "m2_noise2.csv")).
#' Note: the name of each method should be specified. This name will be used for
#' generating the label when plotting.
#' @param condition a vector contains the condition labels corresponding to the predicted
#' proportion file.
#' For example, c("noise0", "noise1", "noise2").
#' @param method One of the c("mae", "rmse", "mape", "kendall", "pearson", "spearman").
#' Note: for mape, cell types with real proportion 0 will be ignored.
#' @param method2 One of the c("mae", "rmse", "mape", "kendall", "pearson", "spearman").
#' Generating the second metric in cheatmap.
#' @param type One of the c("sample", "celltpye", "all").
#' Note: method2 will be disabled when figure is in c("boxplot", "barplot", "heatmap")
#' @param figure One of the c("boxplot", "barplot", "heatmap", "cheatmap")
#' @param errbar error bar type for barplot. One of the c("SD", "SE"), default: SE.
#' SD: Standard Deviation
#' SE: Standard Error
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
#'
#' @export
#'
plot_multiple2 <- function(actual,
                           predicted,
                           condition = NULL,
                           method = NULL,
                           method2 = NULL,
                           type = "sample",
                           figure = "boxplot",
                           errbar = NULL) {

    if (!(method %in% c("mape", "mae", "rmse", "pearson", "spearman", "kendall"))) {
        stop("Parameter method must be one of c('mae', 'mape', 'rmse', 'pearson', 'spearman', 'kendall').")
    }

    if (is.null(condition)) {
        stop("Parameter 'condition' must be provided!")
    }

    # get real proportion
    p_true <- getInput(data = actual, name = "actual")

    if (figure %in% c("boxplot", "barplot")) {
        fi <- list()
        for (name in names(predicted)) {  # name is the method name
            pred_files <- predicted[[name]]
            res <- list()
            for (i in seq(length(pred_files))) {
                tmp_condition <- condition[i]
                p_pred <- getInput(data = pred_files[i], name = "predicted")
                res[[tmp_condition]] <- regMetrics(actual = p_true,
                                                   predicted = p_pred,
                                                   method = method,
                                                   type = type)
            }
            fi[[name]] <- do.call(rbind, res)
        }

    } else if (figure %in% c("heatmap", "cheatmap")) {
        fi <- list()
        for (name in names(predicted)) {  # name is the method name
            pred_files <- predicted[[name]]
            res <- list()
            for (i in seq(length(pred_files))) {
                tmp_condition <- condition[i]
                p_pred <- getInput(data = pred_files[i], name = "predicted")
                res[[tmp_condition]] <- regMetrics(actual = p_true,
                                                   predicted = p_pred,
                                                   method = method,
                                                   type = "all")
            }
            fi[[name]] <- do.call(rbind, res)
        }

        fi2 <- list()
        if (!(is.null(method2))) {
            for (name in names(predicted)) {  # name is the method name
                pred_files <- predicted[[name]]
                res <- list()
                for (i in seq(length(pred_files))) {
                    tmp_condition <- condition[i]
                    p_pred <- getInput(data = pred_files[i], name = "predicted")
                    res[[tmp_condition]] <- regMetrics(actual = p_true,
                                                       predicted = p_pred,
                                                       method = method2,
                                                       type = "all")
                }
                fi2[[name]] <- do.call(rbind, res)
            }
        }

    } else {
        stop(" Parameter figure is not valid! ")
    }

    if (figure == "boxplot") {
        if (type %in% c("sample", "celltype")) {
            title_mess <- paste("Boxplot for", type, method, sep = " ")
        } else {
            stop(" Parameter 'type' must be 'sample' or 'celltype'. ")
        }
        df <- melt(fi)
        colnames(df) <- c("condition", type, "value", "Method")

        p <- ggplot(df, aes(x = .data$condition, y = .data$value, fill = .data$Method)) +
            # stat_boxplot(geom = "errorbar", width = 0.25) +
            geom_boxplot(width=0.4, outlier.shape = 19, outlier.size = 0.3) +
            theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  plot.title = element_text(size = 18, hjust = 0.5),
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5)) +
            labs(title = title_mess, x = "condition", y = method, fill = "Method")

    } else if (figure == "barplot") {
        if (type %in% c("sample", "celltype")) {
            title_mess <- paste("Barplot for", type, method, sep = " ")
        } else {
            stop(" Parameter 'type' must be 'sample' or 'celltype'. ")
        }
        y_mess <- method

        df <- as.data.frame(matrix(ncol = 4, nrow = 0))
        colnames(df) <- c("condition", "Method", "mean", "sd")
        for (this_deMethod in names(fi)) {
            this_data <- fi[[this_deMethod]]
            this_data_mean <- apply(X = this_data, MARGIN = 1, FUN = mean)
            this_data_sd <- apply(X = this_data, MARGIN = 1, FUN = stats::sd)
            this_df <- data.frame(mean = this_data_mean,
                                  sd = this_data_sd,
                                  Method = this_deMethod)
            this_df <- tibble::rownames_to_column(this_df, "condition")
            df <- rbind(df, this_df)
        }

        p <- ggplot(data = df, aes(x = as.factor(.data$condition), y = .data$mean, fill = .data$Method)) +
            geom_bar(position=position_dodge(), stat="identity", colour='black', alpha = 0.7) +
            geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.4, position = position_dodge(0.9), alpha = 0.9) +
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
            ylab(y_mess)
    } else if (figure == "heatmap") {
        writeLines(" Parameter 'figure' is set to 'heatmap', Parameter 'type' will be ignored. ")

        for (name in names(fi)) {
            colnames(fi[[name]]) <- c(name)
        }

        df <- as.matrix(cbind.data.frame(fi))
        df <- round(x = df, digits = 4)
        df <- melt(df)
        colnames(df) <- c("Condition", "Method", "value")

        p <- ggplot(df, aes(x = .data$Condition, y = .data$Method, fill = value)) +
            geom_tile(color = "black") +
            geom_text(aes(label = value), color = "white") +
            coord_fixed() +
            scale_fill_gradient(low = "blue", high = "red") +
            guides(fill = guide_colourbar(title = method))


    } else if (figure == "cheatmap") {
        writeLines(" Parameter 'figure' is set to 'heatmap', Parameter 'type' will be ignored. ")

        for (name in names(fi)) {
            colnames(fi[[name]]) <- c(name)
        }
        df1 <- as.matrix(cbind.data.frame(fi))
        df1 <- round(x = df1, digits = 4)
        df1 <- melt(df1)
        colnames(df1) <- c("Condition", "Method", "value1")


        for (name in names(fi2)) {
            colnames(fi2[[name]]) <- c(name)
        }
        df2 <- as.matrix(cbind.data.frame(fi2))
        df2 <- round(x = df2, digits = 4)
        df2 <- melt(df2)
        colnames(df2) <- c("Condition", "Method", "value2")

        df <- merge(x = df1, y = df2, by = c("Condition", "Method"))

        p <- ggplot(df, aes(x = .data$Condition,
                            y = forcats::fct_rev(.data$Method),
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
            labs(size = method2,
                 fill = method,
                 x = NULL,
                 y = NULL)

    } else {
        stop("Parameter figure is invalid!")
    }

    return(list(data = df, plot = p))

}






















