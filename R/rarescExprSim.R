#' @title Generate in silico mixture expression matrix based on scRNA-seq with rare component.
#'
#' @description Generate in silico mixture expression matrix based on the
#' internal scRNA-seq database. The internal RNA-seq database is collected from
#' PMID:29474909. All the samples are passed the quality filter.
#'
#' @param p_rare A vector of proportions.
#' Default: c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)
#' Note: Every cell type will be treated as rare component.
#' For example, if 8 cell types need to be tested, this function will generate
#' 7 * 8 = 56 samples. 7 means 7 rare proportions and 8 means 8 cell types.
#' @param p Proportion of sample in train set, default: 2 / 3.
#' @param transform "TPM", "CPM" or "NO". Transform the data into TPM, CPM or in raw counts.
#' The suffix will be added to the file name.
#' @param outputPath output file save path.
#' @param bulk_name Pseudo bulk output file name.
#' @param ref_bulk_name reference output file name in csv, for the deconvolution method based on bulk data.
#' @param ref_cell_number number of cells of each cell type for generating reference, default: 1000.
#' @param ref_sc_name File name for all data in train set in csv.
#' This data can be used for differential gene analysis.
#' @param ref_sc_label cell labels for ref_sc_name.
#' @param prop_name simulated proportion file name in csv.
#' @param type 'mouse_tissue' or 'human_PBMC'. Default: 'human_PBMC'
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
#' rarescExprSim()
#'
rarescExprSim <- function (cell_number = 2000,
                           p = 2 / 3,
                           transform = "TPM",
                           outputPath = NULL,
                           bulk_name = "scPBMC_gene_expr.csv",
                           ref_bulk_name = "scPBMC_ref.csv",
                           ref_cell_number = 1000,
                           ref_sc_name = "scPBMC_ref_sc.csv",
                           ref_sc_label = "scPBMC_ref_sc_label.csv",
                           prop_name = "scPBMC_prop.csv",
                           type = 'human_PBMC',
                           p_rare = c(0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05)) {
    writeLines("Loading scRNA-seq data......")

    if(type == 'mouse_tissue') {
        raw_data <- readRDS(file = system.file(package="Deconer", "extdata", "scRNAseq_matrix_mouse.rds"))
        this.celltypes <- c("FetalStomach", "FetalLung", "FetalLiver", "FetalKidney",
                            "FetalIntestine", "FetalBrain", "Female.fetal.Gonad")
    } else if (type == 'human_PBMC') {
        raw_data <- readRDS(file = system.file(package="Deconer", "extdata", "scRNAseq_matrix_PBMC.rds"))
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

    cell_idx_list_train <- c()
    cell_idx_list_test <- c()

    # separate train and test cells
    for (ct in this.celltypes) {
        mess <- paste("Now, split cell type", ct, "into train and test set......", sep = " ")
        writeLines(mess)
        idx <- which(grepl(pattern = ct, x = colnames(data), fixed = TRUE) == TRUE)
        ct_cellnumber <- length(idx)
        ct_trainNum <- as.integer(ct_cellnumber * p)
        ct_testNum <- as.integer(ct_cellnumber - ct_trainNum)

        train_idx <- sample(x = idx, size = ct_trainNum, replace = FALSE)
        test_idx <- setdiff(idx, train_idx)

        cell_idx_list_train[[ct]] <- train_idx
        cell_idx_list_test[[ct]] <- test_idx
    }

    # creating pseudo bulk data
    writeLines("Creating simulated mixture samples......")
    n_sample <- length(p_rare) * length(this.celltypes)
    mess <- paste("There are",
                  length(this.celltypes),
                  "cell types will be tested, and each cell type will generate",
                  length(p_rare), "sample......")
    writeLines(mess)

    prop <- matrix(data = sample(x = 10000, size = n_sample * length(this.celltypes), replace = TRUE),
                   nrow = length(this.celltypes), ncol = n_sample)
    for (ii in seq(length(this.celltypes))) {
        for (jj in seq(length(p_rare))) {
            idx_column <- jj + (ii - 1) * length(p_rare)
            this_sample <- prop[, idx_column]
            other_p <- this_sample[-ii]
            other_p <- (other_p / sum(other_p)) * (1 - p_rare[jj])
            prop[, idx_column] <- append(x = other_p, values = p_rare[jj], (ii - 1))
        }
    }
    colnames(prop) <- paste("S", seq(n_sample), sep = "")
    rownames(prop) <- this.celltypes

    # convert to integer matrix
    cell_number_matrix <- prop * cell_number
    mode(cell_number_matrix) <- "integer"

    # # re-norm prop
    # prop <- apply(X = prop, MARGIN = 2, FUN = Deconer:::v_norm)


    ## generatebulk data
    writeLines("Now, generating pseudo bulk data......")

    pseudo_bulk_expr <- list()

    # creating processing bar
    pb <- progress_bar$new(total = n_sample)

    for(i in colnames(cell_number_matrix)) {
        pb$tick()
        selected_cell_idx <- c()
        this_sample_cell_number <- cell_number_matrix[, i]
        for (ct in this.celltypes) {  # sampling single cells
            this_celltype_cell_number <- this_sample_cell_number[[ct]]
            this_celltype_train_idx <- cell_idx_list_test[[ct]]

            if (length(this_celltype_train_idx) < this_celltype_cell_number) {
                this_celltype_selected_idx <- sample(x = this_celltype_train_idx,
                                                     size = this_celltype_cell_number,
                                                     replace = TRUE)
            } else {
                this_celltype_selected_idx <- sample(x = this_celltype_train_idx,
                                                     size = this_celltype_cell_number,
                                                     replace = FALSE)
            }
            selected_cell_idx <- c(selected_cell_idx, this_celltype_selected_idx)
        }
        pseudo_bulk_expr[[i]] <- apply(X = data[, selected_cell_idx], MARGIN = 1, FUN = sum)

    }

    pseudo_bulk_expr <- as.data.frame(do.call(cbind, pseudo_bulk_expr))


    ## creating reference data
    writeLines("Now, generating reference data......")
    pb <- progress_bar$new(total = length(this.celltypes))
    reference_expr <- list()
    for (ct in this.celltypes) {  # sampling single cells
        pb$tick()
        this_celltype_cell_number <- ref_cell_number
        this_celltype_train_idx <- cell_idx_list_train[[ct]]

        if (length(this_celltype_train_idx) < this_celltype_cell_number) {
            this_celltype_selected_idx <- sample(x = this_celltype_train_idx,
                                                 size = this_celltype_cell_number,
                                                 replace = TRUE)
        } else {
            this_celltype_selected_idx <- sample(x = this_celltype_train_idx,
                                                 size = this_celltype_cell_number,
                                                 replace = FALSE)
        }
        reference_expr[[ct]] <- apply(X = data[, this_celltype_selected_idx], MARGIN = 1, FUN = sum)
    }

    reference_expr <- as.data.frame(do.call(cbind, reference_expr))

    # check the gene order
    if (all(rownames(pseudo_bulk_expr) == rownames(data))) {
        writeLines("pseudo bulk data consistency checking: PASSED!")
    } else {
        stop("pseudo bulk data consistency checking: FAILED!")
    }

    if (all(rownames(reference_expr) == rownames(data))) {
        writeLines("reference data consistency checking: PASSED!")
    } else {
        stop("reference data consistency checking: FAILED!")
    }

    # convert data
    if(transform == "TPM"){
        writeLines("Transform data into TPM......")
        pseudo_bulk_expr_transformed <- Deconer:::TPM(data = Deconer:::merge.all(pseudo_bulk_expr, raw_data["Length"]))
        reference_expr_transformed <- Deconer:::TPM(data = Deconer:::merge.all(reference_expr, raw_data["Length"]))
    }else if(transform == "CPM"){
        writeLines("Transform data into CPM......")
        pseudo_bulk_expr_transformed <- Deconer:::CPM(data = Deconer:::merge.all(pseudo_bulk_expr, raw_data["Length"]))
        reference_expr_transformed <- Deconer:::CPM(data = Deconer:::merge.all(reference_expr, raw_data["Length"]))
    } else {
        writeLines("Note: data transformation is not specified!")
    }

    ## creating reference in single cell format
    reference_celltype_vector <- c()
    reference_idx_vector <- c()
    for (ct in this.celltypes) {  # sampling single cells
        this_celltype_train_idx <- cell_idx_list_train[[ct]]
        reference_celltype_vector <- c(reference_celltype_vector, rep(ct, length(this_celltype_train_idx)))
        reference_idx_vector <- c(reference_idx_vector, this_celltype_train_idx)
    }
    reference_celltype_df <- data.frame(CellType = reference_celltype_vector)

    reference_expr_sc <- data[, reference_idx_vector]

    fwrite(x = as.data.frame(reference_expr),
           file = file.path(outputPath, ref_bulk_name),
           sep = ",",
           row.names = TRUE,
           quote = FALSE)
    fwrite(x = as.data.frame(prop),
           file = file.path(outputPath, prop_name),
           sep = ",",
           row.names = TRUE,
           quote = FALSE)
    fwrite(x = as.data.frame(pseudo_bulk_expr),
           file = file.path(outputPath, bulk_name),
           sep = ",",
           row.names = TRUE,
           quote = FALSE)

    if (transform %in% c("TPM", "CPM")) {
        tmp_name <- file_path_sans_ext(bulk_name)
        bulk_name_transformed <- paste0(tmp_name, "_", transform, ".csv")
        fwrite(x = as.data.frame(pseudo_bulk_expr_transformed),
               file = file.path(outputPath, bulk_name_transformed),
               sep = ",",
               row.names = TRUE,
               quote = FALSE)

        tmp_name <- file_path_sans_ext(ref_bulk_name)
        ref_bulk_name_transformed <- paste0(tmp_name, "_", transform, ".csv")
        fwrite(x = as.data.frame(reference_expr_transformed),
               file = file.path(outputPath, ref_bulk_name_transformed),
               sep = ",",
               row.names = TRUE,
               quote = FALSE)
    }

    if (!is.null(ref_sc_name)) {
        writeLines("Output raw train counts......")
        fwrite(x = as.data.frame(reference_expr_sc),
               file = file.path(outputPath, ref_sc_name),
               sep = ",",
               row.names = TRUE,
               quote = FALSE)
        fwrite(x = reference_celltype_df,
               file = file.path(outputPath, ref_sc_label),
               sep = "\n",
               row.names = FALSE,
               quote = FALSE)
    }
}



