#' @title Generate in silico expression dataset with unknown component.
#'
#' @description Generate in silico mixture expression matrix based on the
#' internal RNA-seq database. This function is different from the function
#' \code{\link{exprSim}}. \code{\link{exprSim}} generates all the proportion
#' randomly, while this function takes one cell type as unknown component.
#' The internal RNA-seq database is collected from multiple studies.
#' All the samples are passed the quality filter. We provide two types of
#' simulation called "coarse" and "fine". This idea is from DREAM Challenge Tumor
#' Deconvolution problem.
#' \url{https://www.synapse.org/#!Synapse:syn15589870/wiki/}.
#'
#' @param unknown a numeric vector defines the proportion of unknown content.
#' Default: c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2)
#' @param n_sample Sample number to be generated, default: 50.
#' @param p Proportion of sample in train set, default: 0.6.
#' @param type "coarse" or "fine".
#' "coarse" means the simulation will be performed in a coarse level.
#' Only 8 cell types will be used, including ("B.cells", "CD4.T.cells",
#' "CD8.T.cells", "endothelial.cells", "macrophages", "monocytes","neutrophils", "NK.cells").
#' "fine" means the simulation will be performed in a fine level.
#' 14 cell types will be used, including ("memory.B.cells", "naive.B.cells", "memory.CD4.T.cells",
#' "naive.CD4.T.cells", "regulatory.T.cells", "memory.CD8.T.cells", "naive.CD8.T.cells", "NK.cells",
#' "neutrophils", "monocytes", "myeloid.dendritic.cells", "macrophages", "fibroblasts", "endothelial.cells")
#' @param transform "TPM", "CPM" or "NO". Transform the data into TPM, CPM or in counts.
#' @param outputPath output file save path.
#' @param mix_name mixture output file name.
#' @param ref_name reference output file name in csv.
#' @param prop_name simulated proportion file name in csv, unknown component
#' will be removed.
#' @param refVar_name reference variance file name in csv.
#' @param train_name file name for all data in train set in csv.
#' This data can be used for differential gene analysis.
#'
#' @return All the information will be written in the output path.
#'
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom utils write.csv
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
#' @examples
#' res <- unExprSim()
#'
unExprSim <- function(unknown = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2),
                      n_sample = 50,
                      p = 0.6,
                      type = "coarse",
                      transform = "TPM",
                      outputPath = NULL,
                      mix_name = "coarse_gene_expr.csv",
                      ref_name = "coarse_ref.csv",
                      prop_name = "coarse_prop.csv",
                      refVar_name = NULL,
                      train_name = NULL) {
    if (is.null(unknown)) {
        stop( "Parameter unknown must be provided!" )
    }

    writeLines("Loading RNA-seq data......")
    attr <- readRDS(file = system.file(package="decone", "extdata", "RNAseq_attr.rds"))
    raw_data <- readRDS(file = system.file(package="decone", "extdata", "RNAseq_matrix.rds"))

    if (is.null(outputPath)) {
        writeLines("output path is not specified, all the file will be saved in work directory.")
        outputPath <- "."
    } else {
        if (!dir.exists(outputPath)) {
            dir.create(outputPath)
        }
    }

    # data_counts is not transformed data
    data <- raw_data[ , -which(names(raw_data) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length"))]
    data_counts <- data

    if(transform == "TPM"){
        writeLines("Transform data into TPM......")
        data <- decone:::TPM(data = decone:::merge.all(data, raw_data["Length"]))
    }else if(transform == "CPM"){
        writeLines("Transform data into CPM......")
        data <- decone:::CPM(data = merge.all(decone:::data, raw_data["Length"]))
    } else {
        writeLines("Note: data transformation is not specified!")
    }

    if (type == "coarse") {
        this.celltypes <- c("B.cells", "CD4.T.cells", "CD8.T.cells", "endothelial.cells",
                            "macrophages", "monocytes","neutrophils", "NK.cells")

    } else if (type == "fine") {
        this.celltypes <- c("memory.B.cells", "naive.B.cells", "memory.CD4.T.cells", "naive.CD4.T.cells",
                            "regulatory.T.cells", "memory.CD8.T.cells", "naive.CD8.T.cells", "NK.cells",
                            "neutrophils", "monocytes", "myeloid.dendritic.cells", "macrophages",
                            "fibroblasts", "endothelial.cells")
    } else {
        stop("Parameter type must be coarse or fine.")
    }

    writeLines("Performing train test split......")
    this.attr <- attr[which(attr$Required_Cell_Type %in% this.celltypes), ]
    this.list <- decone:::multiClassSample(df = this.attr, p = p)

    this.train <- this.list[["train"]]
    this.test <- this.list[["test"]]

    writeLines("Creating train set......")
    this.train_matrix <- data %>%
        subset(subset = TRUE, select = this.list[["train"]])

    writeLines("Creating test set......")
    this.test_matrix <- data %>%
        subset(subset = TRUE, select = this.list[["test"]])


    writeLines("Generating train reference data......")
    this.train <- decone:::refGenerator(df = this.train_matrix,
                                        attr = this.attr,
                                        cellTypes = this.celltypes)
    thisref.train <- this.train$ref
    thisstd.train <- this.train$std

    writeLines("Generating test reference data......")
    this.test <- decone:::refGenerator(df = this.test_matrix,
                                       attr = this.attr,
                                       cellTypes = this.celltypes)
    thisref.test <- this.test$ref
    thisstd.test <- this.test$std


    writeLines("Creating simulated mixture samples......")

    writeLines("Now, generating samples without unknown content......")
    prop <- matrix(data = sample(x = 1000, size = n_sample * length(this.celltypes), replace = TRUE),
                   nrow = length(this.celltypes), ncol = n_sample)
    prop <- apply(X = prop, MARGIN = 2, FUN = v_norm)
    colnames(prop) <- paste("S", seq(n_sample), sep = "")
    rownames(prop) <- colnames(thisref.test)
    this.mix <- thisref.test %*% prop

    write.csv(x = as.data.frame(thisref.train),
              file = file.path(outputPath, ref_name),
              row.names = TRUE)
    write.csv(x = prop,
              file = file.path(outputPath, prop_name),
              row.names = TRUE)
    write.csv(x = this.mix,
              file = file.path(outputPath, mix_name),
              row.names = TRUE)

    if (!is.null(train_name)) {
        writeLines("Output all samples in train set......")
        this.train_counts <- data_counts %>%
            subset(subset = TRUE, select = this.list[["train"]])
        write.csv(x = this.train_counts,
                  file = file.path(outputPath, train_name),
                  row.names = TRUE)
    }

    if (!is.null(refVar_name)) {
        writeLines("Output reference variance file......")
        write.csv(x = thisstd.train,
                  file = file.path(outputPath, refVar_name),
                  row.names = TRUE)
    }



    for (un_p in unknown) {
        mess <- paste("Now, generating samples with", un_p, "unknown content......", sep = " ")
        writeLines(mess)
        prop <- matrix(data = sample(x = 1000, size = n_sample * length(this.celltypes), replace = TRUE),
                       nrow = length(this.celltypes), ncol = n_sample)
        prop <- apply(X = prop, MARGIN = 2, FUN = v_norm, scale = (1 - un_p))
        colnames(prop) <- paste("S", seq(n_sample), sep = "")
        rownames(prop) <- colnames(thisref.test)
        this.mix <- thisref.test %*% prop

        tmp_name <- paste0(file_path_sans_ext(mix_name), "_un_", un_p, ".csv")
        write.csv(x = this.mix,
                  file = file.path(outputPath, tmp_name),
                  row.names = TRUE)
    }
}
