#' @title Boxplot of deconvolution results for multiple test.
#'
#' @description Boxplot for illustrating the deconvolution performance with the
#' noised in silico data.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a vector contains all the files with the predicted proportions
#' in different noise level. Must be in csv file with row and column names.
#' row: cell types, column: samples.
#' For example, c("noise0.csv", "noise1.csv", "noise2.csv").
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: base file name in parameter 'predicted'.
#' For example, c("noise0", "noise1", "noise2").
#' @param method "rmse", "mape", "mae", "pearson" or "spearman".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#' @param title Boxplot title in character.
#'
#' @return The computed metrics and plot data.
#'
#' @importFrom Metrics mape rmse mae
#' @importFrom ggplot2 ggplot geom_boxplot aes stat_boxplot theme
#' ggtitle xlab ylab element_text element_rect element_blank .data
#' @importFrom tools file_path_sans_ext
#' @importFrom reshape2 melt
#' @importFrom utils read.csv
#'
#' @export
#'
#' @examples
#' res <- pseudoData(type = 2)
#' res <-boxplot_NGrad(actual = res$actual, predicted = res$predicted,
#' label = res$noise_level, method = "rmse")
#'
boxplot_NGrad <- function(actual,
                              predicted,
                              label = NULL,
                              method,
                              title = "Boxplot") {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (!(method %in% c("mape", "mae", "rmse", "pearson", "spearman"))) {
        stop("Parameter must be mape, mae, rmse, pearson or spearman.")
    }

    if (is.null(label)) {
        label <- file_path_sans_ext(basename(predicted))
    }

    if (!(length(label) == length(predicted))) {
        stop("Different length between parameter predicted and label.")
    }

    res <- list()

    for (i in seq(length(predicted))) {
        this.label <- label[i]
        p_pred <- as.matrix(read.csv(file = predicted[i],
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))

        res[[this.label]] <- regMetrics(actual = p_true,
                                        predicted = p_pred,
                                        method = method)
    }

    df <- do.call(cbind, res)
    df1 <- melt(data = df)

    p <- ggplot(data = df1, aes(x = .data$Var2, y = .data$value)) +
        # stat_boxplot(geom = "errorbar", width = 0.25) +
        geom_boxplot(width=0.4, outlier.shape = 19) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, size = 1),
              plot.title = element_text(size = 18, hjust = 0.5),
              axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
              axis.title.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
              axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5)) +
        ggtitle(title) +
        xlab("Noise Level") +
        ylab(method)

    return(list(data = df1, plot = p))

}



#' @title Comparing boxplot for noise testing.
#'
#' @description Boxplot for illustrating the deconvlition performance with the
#' noised in silico data. this function is used for comparing multiple methods.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a list of vectors. Each vector contains all the files with
#' the predicted proportions in different noise level. Must be in csv file with
#' row and column names.
#' row: cell types, column: samples.
#' For example, list(method1 = c("m1_noise0.csv", "m1_noise1.csv", "m1_noise2.csv"),
#' method2 = c("m2_noise0.csv", "m2_noise1.csv", "m2_noise2.csv")).
#' Note: names for each method will be used for assign colors.
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: c("condition1", "condition2", "condition3", ...)
#' @param method "rmse", "mape", "mae", "pearson" or "spearman".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#' @param title Boxplot title in character. Default: "Boxplot".
#'
#' @return The computed metrics and plot data.
#'
#' @importFrom Metrics mape rmse mae
#' @importFrom ggplot2 ggplot geom_boxplot aes stat_boxplot theme
#' ggtitle xlab ylab element_text element_rect element_blank labs .data
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' res <- pseudoData(type = 3, outputPath = "test_file")
#' res <- boxplot_NcrossCompare(actual = res$actual, predicted = res$predicted,
#' label = res$noise_level, method = "pearson", title = "pearson")
#'
boxplot_NcrossCompare <- function(actual,
                                 predicted,
                                 label = NULL,
                                 method,
                                 title = "Boxplot") {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (!(method %in% c("mape", "mae", "rmse", "pearson", "spearman"))) {
        stop("Parameter must be mape, mae, rmse, pearson, spearman.")
    }

    if (is.null(label)) {
        label <- paste0("condition", seq(length(predicted[[1]])))
    }

    fi <- list()

    for (tmp_name in names(predicted)) {
        pred_files <- predicted[[tmp_name]]

        if (!(length(label) == length(pred_files))) {
            stop("File number is inconsistent in parameter 'predicted'.")
        }

        res <- list()
        for (i in seq(length(pred_files))) {
            tmp_label <- label[i]
            p_pred <- as.matrix(read.csv(file = pred_files[i],
                                         header = TRUE,
                                         row.names = 1,
                                         encoding = "UTF-8"))
            res[[tmp_label]] <- regMetrics(actual = p_true, predicted = p_pred, method = method)
        }
        fi[[tmp_name]] <- do.call(cbind, res)
    }

    fi <- melt(fi)

    p <- ggplot(fi, aes(x = .data$Var2, y = .data$value, fill = .data$L1)) +
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
        labs(title = title, x = "Noise Level", y = method, fill = "Method")

    return(list(data = fi, plot = p))
}



#' @title Plot heatmap for deconvolution results.
#'
#' @description Generate heatmap comparison for different methods.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a list of vectors. Each vector contains all the files with
#' the predicted proportions in different noise level. Must be in csv file with
#' row and column names.
#' row: cell types, column: samples.
#' For example, list(method1 = c("m1_noise0.csv", "m1_noise1.csv", "m1_noise2.csv"),
#' method2 = c("m2_noise0.csv", "m2_noise1.csv", "m2_noise2.csv")).
#' Note: names for each method will be used for assign colors.
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: c("condition1", "condition2", "condition3", ...)
#' @param method "rmse", "mape", "mae", "pearson" or "spearman".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#'
#' @return The computed metrics and plot data.
#'
#' @importFrom Metrics mape rmse mae
#' @importFrom stats cor
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices heat.colors
#'
#' @export
#'
#' @examples
#' res <- pseudoData(type = 3, outputPath = "test_file")
#' res <- heatmap_NcrossCompare(actual = res$actual, predicted = res$predicted,
#' label = res$noise_level, method = "mape")
#'
heatmap_NcrossCompare <- function(actual,
                                 predicted,
                                 label = NULL,
                                 method) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (!(method %in% c("mape", "mae", "rmse", "pearson", "spearman"))) {
        stop("Parameter must be mape, mae, rmse, pearson or spearman.")
    }

    if (is.null(label)) {
        label <- paste0("condition", seq(length(predicted[[1]])))
    }

    fi <- list()

    for (name in names(predicted)) {
        pred_files <- predicted[[name]]
        res <- list()
        for (i in seq(length(pred_files))) {
            tmp_label <- label[i]
            p_pred <- as.matrix(read.csv(file = pred_files[i],
                                         header = TRUE, row.names = 1,
                                         encoding = "UTF-8"))

            if (method %in% c("mape", "mae", "rmse")) {
                res[[tmp_label]] <- mean(regMetrics(actual = p_true,
                                                    predicted = p_pred,
                                                    method = method))
            } else {
                method <- paste0(method, "2")
                res[[tmp_label]] <- regMetrics(actual = p_true,
                                               predicted = p_pred,
                                               method = method)
            }
        }
        fi[[name]] <- unlist(res)
    }

    fi <- do.call(rbind, fi)

    p <- pheatmap(mat = fi,
                  display_numbers = TRUE,
                  color = rev(heat.colors(n = ncol(fi))),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  cellwidth = 40,
                  cellheight = 20,
                  fontsize = 16)

    return(list(data = fi, plot = p))
}



#' @title Plot heatmap for different cell types.
#'
#' @description Generate heatmap comparison for different cell types.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a vector contains all the files with the predicted proportions
#' in different noise level. Must be in csv file with row and column names.
#' row: cell types, column: samples.
#' For example, c("noise0.csv", "noise1.csv", "noise2.csv").
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: c("condition1", "condition2", "condition3", ...)
#' @param method "rmse", "mape", "mae", "pearson" or "spearman".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#'
#' @return The computed metrics and plot data.
#'
#' @importFrom pheatmap pheatmap
#'
#' @export
#'
#' @examples
#' res <- pseudoData(type = 2)
#' res <- heatmap_NGradCT(actual = res$actual, predicted = res$predicted,
#' label = res$noise_level, method = "rmse")
#'
heatmap_NGradCT <- function (actual,
                             predicted,
                             label = NULL,
                             method) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (!(method %in% c("mape", "mae", "rmse", "pearson", "spearman"))) {
        stop("Parameter must be mape, mae, rmse, pearson or spearman.")
    }

    if (is.null(label)) {
        label <- paste0("condition", seq(length(predicted[[1]])))
    }

    res <- list()
    for (i in seq(length(predicted))) {
        this.label <- label[i]
        p_pred <- as.matrix(read.csv(file = predicted[i],
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))

        res[[this.label]] <- regMetrics(actual = t(p_true),
                                        predicted = t(p_pred),
                                        method = method)
    }

    df <- do.call(cbind, res)

    p <- pheatmap(mat = df,
                  display_numbers = TRUE,
                  color = rev(heat.colors(n = ncol(df))),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  cellwidth = 40,
                  cellheight = 20,
                  fontsize = 16)

    return(list(data = df, plot = p))
}




#' @title Plot circle heatmap for deconvolution results.
#'
#' @description Usually, only one metric cannot reveal the deconvolution
#' performance well. Therefore, the circle heatmap which contains two different
#' dimension information can be used for a better illustration.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a list of vectors. Each vector contains all the files with
#' the predicted proportions in different noise level. Must be in csv file with
#' row and column names.
#' row: cell types, column: samples.
#' For example, list(method1 = c("m1_noise0.csv", "m1_noise1.csv", "m1_noise2.csv"),
#' method2 = c("m2_noise0.csv", "m2_noise1.csv", "m2_noise2.csv")).
#' Note: names for each method will be used for assign colors.
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: c("condition1", "condition2", "condition3", ...)
#' @param method1 "rmse", "mape", "mae", "pearson" or "spearman".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#' @param method2 "rmse", "mape", "mae", "pearson" or "spearman".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#'
#' @return The computed metrics and plot data.
#'
#' @importFrom reshape2 melt
#' @importFrom forcats fct_rev
#' @importFrom ggplot2 ggplot geom_point scale_x_discrete scale_radius
#' scale_fill_gradient theme_minimal theme element_text element_blank
#' guides guide_legend guide_colorbar .data
#' @importFrom stringr str_remove
#'
#' @export
#'
#' @examples
#' res <- pseudoData(type = 3, outputPath = "test_file")
#' res <- cheatmap_NcrossCompare(actual = res$actual, predicted = res$predicted,
#' label = res$noise_level,  method1 = "pearson", method2 = "rmse")
#'
cheatmap_NcrossCompare <- function(actual, predicted, label = NULL, method1, method2) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (!(method1 %in% c("mape", "mae", "rmse", "pearson", "spearman"))) {
        stop("Parameter must be mape, mae, rmse, pearson or spearman.")
    }

    if (!(method2 %in% c("mape", "mae", "rmse", "pearson", "spearman"))) {
        stop("Parameter must be mape, mae, rmse, pearson or spearman.")
    }

    if (method1 %in% c("pearson", "spearman")) {
        method1 <- paste0(method1, "2")
    }

    if (method2 %in% c("pearson", "spearman")) {
        method2 <- paste0(method2, "2")
    }

    if (is.null(label)) {
        label <- paste0("condition", seq(length(predicted[[1]])))
    }

    fi1 <- list()
    fi2 <- list()
    for (name in names(predicted)) {
        pred_files <- predicted[[name]]
        res1 <- list()
        res2 <- list()
        for (i in seq(length(pred_files))) {
            this_label <- label[i]
            p_pred <- as.matrix(read.csv(file = pred_files[i],
                                         header = TRUE, row.names = 1,
                                         encoding = "UTF-8"))

            if (method1 %in% c("mape", "mae", "rmse")) {
                res1[[this_label]] <- mean(regMetrics(actual = p_true,
                                                predicted = p_pred,
                                                method = method1))
            } else {
                res1[[this_label]] <- regMetrics(actual = p_true,
                                           predicted = p_pred,
                                           method = method1)
            }

            if (method2 %in% c("mape", "mae", "rmse")) {
                res2[[this_label]] <- mean(regMetrics(actual = p_true,
                                                predicted = p_pred,
                                                method = method2))
            } else {
                res2[[this_label]] <- regMetrics(actual = p_true,
                                           predicted = p_pred,
                                           method = method2)
            }

        }
        fi1[[name]] <- unlist(res1)
        fi2[[name]] <- unlist(res2)
    }
    fi1 <- do.call(rbind, fi1)
    fi2 <- do.call(rbind, fi2)

    a1 <- melt(fi1)
    a2 <- melt(fi2)

    data <- merge(x = a1, y = a2, by = c("Var1", "Var2"))

    if (method1 %in% c("pearson2", "spearman2")) {
        method1 <- str_remove(string = method1, pattern = "2")
        trend_label <- c("Bad", "Good")
    } else {
        trend_label <- c("Good", "Bad")
    }

    if (method2 %in% c("pearson2", "spearman2")) {
        method2 <- str_remove(string = method2, pattern = "2")
    }


    # fill is fi1, size is fi2
    p <- ggplot(data, aes(x = .data$Var2,
                          y = forcats::fct_rev(.data$Var1),
                          fill = .data$value.x,
                          size = .data$value.y)) +
        geom_point(shape = 21,
                   stroke = 0) +
        scale_x_discrete(position = "bottom") +
        scale_radius(range = c(5, 15)) +
        scale_fill_gradient(low = "orange",
                            high = "blue",
                            breaks = c(min(data$value.x), max(data$value.x)),
                            # labels = trend_label,
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
             fill = method1,
             x = NULL,
             y = NULL)

    return(list(data = data, plot = p))
}






