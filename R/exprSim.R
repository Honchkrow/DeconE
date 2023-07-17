#' @title Generate in silico mixture expression matrix
#'
#' @description Generate in silico mixture expression matrix based on the
#' internal RNA-seq database. The internal RNA-seq database is collected from
#' multiple studies. All the samples are passed the quality filter. This
#' function is different from \code{\link{pseudoExpr}}. Counts in
#' \code{\link{pseudoExpr}} is randomly generated from uniform distribution.
#' This function use the real cell type specific expression data to generate
#' mixture data. We provide two types of simulation called "coarse" and "fine".
#' This idea is from DREAM Challenge Tumor Deconvolution problem.
#' \url{https://www.synapse.org/#!Synapse:syn15589870/wiki/}.
#'
#' @param n_sample Sample number to be generated, default: 50.
#' @param p Proportion of sample in train set, default: 2/3.
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
#' @param prop_name simulated proportion file name in csv.
#' @param refVar_name reference variance file name in csv.
#' @param train_name file name for all data in train set in csv.
#' This data can be used for differential gene analysis.
#'
#' @return All the information will be written in the output path.
#'
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom utils write.csv
#'
#' @export
#'
#' @examples
#' exprSim()
#'
exprSim <- function(n_sample = 50,
                    p = 2 / 3,
                    type = "coarse",
                    transform = "TPM",
                    outputPath = NULL,
                    mix_name = "coarse_gene_expr.csv",
                    ref_name = "coarse_ref.csv",
                    prop_name = "coarse_prop.csv",
                    refVar_name = NULL,
                    train_name = NULL) {
    writeLines("Loading RNA-seq data......")
    attr <- readRDS(file = system.file(package="Deconer", "extdata", "RNAseq_attr.rds"))
    raw_data <- readRDS(file = system.file(package="Deconer", "extdata", "RNAseq_matrix.rds"))

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
        data <- TPM(data = merge.all(data, raw_data["Length"]))
    }else if(transform == "CPM"){
        writeLines("Transform data into CPM......")
        data <- CPM(data = merge.all(data, raw_data["Length"]))
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
    this.list <- multiClassSample(df = this.attr, p = p)

    this.train <- this.list[["train"]]
    this.test <- this.list[["test"]]

    writeLines("Creating train set......")
    this.train_matrix <- data %>%
        subset(subset = TRUE, select = this.list[["train"]])

    writeLines("Creating test set......")
    this.test_matrix <- data %>%
        subset(subset = TRUE, select = this.list[["test"]])

    if (!is.null(train_name)) {
        writeLines("Output all samples in train set......")
        this.train_counts <- data_counts %>%
            subset(subset = TRUE, select = this.list[["train"]])
        write.csv(x = this.train_counts,
                  file = file.path(outputPath, train_name),
                  row.names = TRUE)
    }

    writeLines("Generating train reference data......")
    this.train <- refGenerator(df = this.train_matrix,
                               attr = this.attr,
                               cellTypes = this.celltypes)
    thisref.train <- this.train$ref
    thisstd.train <- this.train$std

    writeLines("Generating test reference data......")
    this.test <- refGenerator(df = this.test_matrix,
                              attr = this.attr,
                              cellTypes = this.celltypes)
    thisref.test <- this.test$ref
    thisstd.test <- this.test$std

    if (!is.null(refVar_name)) {
        writeLines("Output reference variance file......")
        write.csv(x = thisstd.train,
                  file = file.path(outputPath, refVar_name),
                  row.names = TRUE)
    }


    writeLines("Creating simulated mixture samples......")
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

}






