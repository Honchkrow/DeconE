#' @title Generate in silico mixture expression matrix based on scRNA-seq
#'
#' @description Generate in silico mixture expression matrix based on the
#' internal scRNA-seq database. The internal RNA-seq database is collected from
#' PMID:29474909. All the samples are passed the quality filter.
#' This function use the real cell type specific expression data to generate
#' mixture data. The cell type include c("FetalStomach", "FetalLung",
#' "FetalLiver", "FetalKidney", "FetalIntestine",
#' "FetalBrain", "Female.fetal.Gonad")
#'
#' @param n_sample Sample number to be generated, default: 50.
#' @param p Proportion of sample in train set, default: 2 / 3.
#' @param transform "TPM", "CPM" or "NO". Transform the data into TPM, CPM or in counts.
#' @param outputPath output file save path.
#' @param mix_name mixture output file name.
#' @param ref_name reference output file name in csv.
#' @param prop_name simulated proportion file name in csv.
#' @param train_name file name for all data in train set in csv.
#' This data can be used for differential gene analysis.
#' @param type 'mouse_tissue' or 'human_PBMC'. Default: 'mouse_tissue'
#'
#' @return All the information will be written in the output path.
#'
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom data.table fwrite
#' @importFrom stringr str_detect
#'
#' @export
#'
#' @examples
#' res <- scExprSim()
#'
scExprSim <- function (n_sample = 50,
                       p = 2 / 3,
                       transform = "TPM",
                       outputPath = NULL,
                       mix_name = "scMouse_gene_expr.csv",
                       ref_name = "scMouse_ref.csv",
                       prop_name = "scMouse_prop.csv",
                       train_name = "scMouse_ref_rawCount.csv",
                       type = 'mouse_tissue') {
    writeLines("Loading scRNA-seq data......")

    if(type == 'mouse_tissue') {
        raw_data <- readRDS(file = system.file(package="decone", "extdata", "scRNAseq_matrix_mouse.rds"))
        this.celltypes <- c("FetalStomach", "FetalLung", "FetalLiver", "FetalKidney",
                            "FetalIntestine", "FetalBrain", "Female.fetal.Gonad")
    } else if (type == 'human_PBMC') {
        raw_data <- readRDS(file = system.file(package="decone", "extdata", "scRNAseq_matrix_PBMC.rds"))
        this.celltypes <- c("intermediate_mono", "CD8+_naïve_T", "mDC", "CD4+_naïve_T", "NK",
                            "memory_B", "CD4+_memory_T", "CD16_mono", "pDC", "naïve_B",
                            "CD8+_activated_T", "CD14_mono", "MAIT")
    } else {
        stop(" Parameter 'type' is not valid! ")
    }

    if (is.null(outputPath)) {
        writeLines("output path is not specified, all the file will be saved in work directory.")
        outputPath <- "."
    } else {
        if (!dir.exists(outputPath)) {
            dir.create(outputPath)
        }
    }

    data <- raw_data[ , -which(colnames(raw_data) %in% c("Length"))]

    train_v <- c()
    test_v <- c()

    train_df <- list()
    test_df <- list()

    for (ct in this.celltypes) {
        mess <- paste("Now, sampling cell type", ct, "......", sep = " ")
        writeLines(mess)
        idx <- which(grepl(pattern = ct, x = colnames(data), fixed = TRUE) == TRUE)
        total_num <- length(idx)
        ct_trainNum <- as.integer(total_num * p)
        ct_testNum <- as.integer(total_num - ct_trainNum)
        train_idx <- sample(x = idx, size = ct_trainNum)
        test_idx <- setdiff(idx, train_idx)

        train_v <- c(train_v, train_idx)
        test_v <- c(test_v, test_idx)

        train_df[[ct]] <- apply(X = data[, train_idx], MARGIN = 1, FUN = sum)
        test_df[[ct]] <- apply(X = data[, test_idx], MARGIN = 1, FUN = sum)
    }

    train_df <- as.data.frame(do.call(cbind, train_df))
    test_df <- as.data.frame(do.call(cbind, test_df))

    if (all(rownames(train_df) == rownames(test_df))) {
        writeLines("Train and test data consistency checking: PASSED!")
    } else {
        stop("Train and test data consistency checking: FAILED!")
    }

    # output for DE
    train_raw <- data[, train_v]
    test_raw <- data[, test_v]

    if(transform == "TPM"){
        writeLines("Transform data into TPM......")
        train_df <- decone:::TPM(data = decone:::merge.all(train_df, raw_data["Length"]))
        test_df <- decone:::TPM(data = decone:::merge.all(test_df, raw_data["Length"]))
    }else if(transform == "CPM"){
        writeLines("Transform data into CPM......")
        train_df <- decone:::CPM(data = decone:::merge.all(train_df, raw_data["Length"]))
        test_df <- decone:::CPM(data = decone:::merge.all(test_df, raw_data["Length"]))
    } else {
        writeLines("Note: data transformation is not specified!")
    }

    writeLines("Creating simulated mixture samples......")
    prop <- matrix(data = sample(x = 10000, size = n_sample * length(this.celltypes), replace = TRUE),
                   nrow = length(this.celltypes), ncol = n_sample)
    prop <- apply(X = prop, MARGIN = 2, FUN = v_norm)
    colnames(prop) <- paste("S", seq(n_sample), sep = "")
    rownames(prop) <- colnames(test_df)
    mix <- test_df %*% prop

    fwrite(x = as.data.frame(train_df),
           file = file.path(outputPath, ref_name),
           sep = ",",
           row.names = TRUE,
           quote = FALSE)
    fwrite(x = as.data.frame(prop),
           file = file.path(outputPath, prop_name),
           sep = ",",
           row.names = TRUE,
           quote = FALSE)
    fwrite(x = as.data.frame(mix),
           file = file.path(outputPath, mix_name),
           sep = ",",
           row.names = TRUE,
           quote = FALSE)

    if (!is.null(train_name)) {
        writeLines("Output raw train counts......")
        fwrite(x = train_raw,
               file = file.path(outputPath, train_name),
               sep = ",",
               row.names = TRUE,
               quote = FALSE)
    }

}



