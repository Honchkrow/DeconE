#' @title Scatter plot for rare component.
#'
#' @description Generate scatter plot for rare component.
#'
#' @param actual The groundtruth proportion of cell types in matrix or csv file.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted The predicted proportion of cell types in matrix or csv file.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param p_rare A vector of proportions. should be the same with function \code{\link{rareExprSim}}.
#' Default: c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)
#' @param celltype TRUE or FALSE. Assign different type and color for different
#' cell types. Default: TRUE.
#'
#' @return The plot data.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_point scale_shape_manual geom_smooth
#' theme element_blank element_rect element_text labs .data
#' @importFrom ggpubr stat_cor
#' @importFrom stats lm
#'
#' @export
#'
#' @examples
#' rareExprSim()
#' # after the simulation, just pass the predicted file
#' # to the parameter "actual" and "predicted"
#'
scatter_R <- function(actual,
                      predicted,
                      p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05),
                      celltype = TRUE) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    } else {
        p_true <- actual
    }

    if (typeof(predicted) == "character") {
        p_pred <- as.matrix(read.csv(file = predicted,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    } else {
        p_pred <- predicted
    }

    if (!all(dim(p_pred) == dim(p_true))) {
        stop("The dimension of input actual and predicted is not consistent.")
    }

    p_len <- length(p_rare)

    pt_list <- list()
    pp_list <- list()

    for (i in seq(nrow(p_true))) {
        ct <- rownames(p_true)[i]
        idx_start <- (i - 1) * p_len + 1
        idx_end <- i * p_len
        print(idx_start)
        print(idx_end)
        pt_list[[ct]] <- p_true[ct, idx_start:idx_end]
        pp_list[[ct]] <- p_pred[ct, idx_start:idx_end]
    }

    pt_m <- do.call(rbind, pt_list)
    colnames(pt_m) <- p_rare
    pp_m <- do.call(rbind, pp_list)
    colnames(pp_m) <- p_rare

    data <- melt(pp_m)

    if (celltype) {
        p <- ggplot(data, aes(x = .data$Var2, y = .data$value)) +
            geom_point(aes(shape = .data$Var1, color = .data$Var1)) +
            scale_shape_manual(values = seq(length(unique(data$Var1)))) +
            geom_smooth(method = lm , color = "red", fill = "#69b3a2", se = TRUE) +
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
        p <- ggplot(data, aes(x = .data$Var2, y = .data$value)) +
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

    return(list(data = data, plot = p))

}


#' @title Plot heatmap for deconvolution results.
#'
#' @description Generate heatmap comparison for different methods.
#'
#' @param actual The groundtruth proportion of cell types in matrix.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a vector contains all the files with the predicted proportions
#' in different method. Must be in csv file with row and column names.
#' row: cell types, column: samples.
#' For example, c("method1.csv", "method2.csv", "method3.csv").
#' @param p_rare A vector of proportions. should be the same with function \code{\link{rareExprSim}}.
#' Default: c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: c("condition1", "condition2", "condition3", ...)
#' @param method "rmse", "mape", "mae".
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
#' rareExprSim()
#' # after the simulation, just pass the predicted file
#' # to the parameter "actual" and "predicted"
#'
heatmap_RcrossCompare <- function (actual,
                                   predicted,
                                   p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05),
                                   label = NULL,
                                   method) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    if (!(method %in% c("mape", "mae", "rmse"))) {
        stop("Parameter must be mape, mae, rmse.")
    }

    if (is.null(label)) {
        label <- paste0("condition", seq(length(predicted[[1]])))
    }

    p_len <- length(p_rare)

    if (! (length(predicted) == length(label))) {
        stop("Parameter predicted and label are not consistant.")
    }

    res <- list()

    for (i in seq(length(predicted))) {

        p_pred <- as.matrix(read.csv(file = predicted[i],
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
        tmp_label <- label[i]

        if (!all(dim(p_pred) == dim(p_true))) {
            mess <- paste("The dimension of input file",
                          predicted[1],
                          "is not consistent.")
            stop(mess)
        }

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

        res[[tmp_label]] <- regMetrics(actual = pt_m,
                                       predicted = pp_m,
                                       method = method)
    }

    df <- do.call(rbind, res)

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





#' @title Plot circle heatmap for deconvolution results rare cell types.
#'
#' @description Generate circle heatmap for different rare cell type proportions.
#' The dataset should be generated from the function \code{\link{rareExprSim}}.
#'
#' @param actual The groundtruth proportion of cell types in matrix or csv file.
#' row: cell types, column: samples. Also can be a csv matrix with row and
#' column names.
#' @param predicted a vector contains all the files with the predicted proportions
#' in different method. Must be in csv file with row and column names.
#' row: cell types, column: samples.
#' For example, c("method1.csv", "method2.csv", "method3.csv").
#' @param p_rare A vector of proportions. should be the same with function \code{\link{rareExprSim}}.
#' Default: c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)
#' @param label a vector contains the label corresponding to the predicted
#' proportion file. Default: c("condition1", "condition2", "condition3", ...)
#' @param method1 "rmse", "mape", "mae".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#' @param method2 "rmse", "mape", "mae".
#' Note: for mape, cell types with read proportion 0 will be ignored.
#'
#' @return The plot data.
#'
#' @importFrom reshape2 melt
#' @importFrom forcats fct_rev
#' @importFrom ggplot2 ggplot geom_point scale_x_discrete scale_radius
#' scale_fill_gradient theme_minimal theme element_text element_blank
#' guides guide_legend guide_colorbar .data
#'
#' @export
#'
#' @examples
#' rareExprSim()
#' # after the simulation, just pass the predicted file
#' # to the parameter "actual" and "predicted"
#'
cheatmap_RcrossCompare <- function (actual,
                                    predicted,
                                    p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05),
                                    label = NULL,
                                    method1,
                                    method2) {
    if (typeof(actual) == "character") {
        p_true <- as.matrix(read.csv(file = actual,
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
    }

    p_len <- length(p_rare)

    if (!(method1 %in% c("mape", "mae", "rmse"))) {
        stop("Parameter must be mape, mae, rmse.")
    }

    if (!(method2 %in% c("mape", "mae", "rmse"))) {
        stop("Parameter must be mape, mae, rmse.")
    }

    if (is.null(label)) {
        label <- paste0("condition", seq(length(predicted[[1]])))
    }


    if (! (length(predicted) == length(label))) {
        stop("Parameter predicted and label are not consistant.")
    }

    fi1 <- list()
    fi2 <- list()
    for (i in seq(length(predicted))) {

        p_pred <- as.matrix(read.csv(file = predicted[i],
                                     header = TRUE,
                                     row.names = 1,
                                     encoding = "UTF-8"))
        tmp_label <- label[i]

        if (!all(dim(p_pred) == dim(p_true))) {
            mess <- paste("The dimension of input file",
                          predicted[1],
                          "is not consistent.")
            stop(mess)
        }

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

        fi1[[tmp_label]] <- regMetrics(actual = pt_m,
                                       predicted = pp_m,
                                       method = method1)
        fi2[[tmp_label]] <- regMetrics(actual = pt_m,
                                       predicted = pp_m,
                                       method = method2)

    }

    fi1 <- do.call(rbind, fi1)
    fi2 <- do.call(rbind, fi2)

    a1 <- melt(fi1)
    a2 <- melt(fi2)

    data <- merge(x = a1, y = a2, by = c("Var1", "Var2"))
    data$Var2 <- as.character(data$Var2)

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








