#' @title Add noise to the simulated expression matrix
#'
#' @description Add different level noise to the simulated expression matrix
#' based on the negative binomial distribution. The related model can be found
#' in this paper (Jin, H., Liu, Z. A benchmark for RNA-seq deconvolution analysis
#'  under dynamic testing environments. Genome Biol 22, 102 (2021).
#'  https://doi.org/10.1186/s13059-021-02290-6).
#'
#' @param exprFile The input expression file. Must be in csv format. Each row is
#' a gene, each column is a sample. rownames and colnames are required. Please check
#' the output of the function \code{\link{exprSim}}.
#' @param outputPath Output path, create if not exists. Default: a new folder
#' based on the exprFile. For example, if the exprFile is "/data/test.txt",
#' the outputPath will be a new folder "/data/test".
#' @param prefix The prefix of the output file.
#' Note: the suffix will be added automatically based on the noise level.
#' @param Pt Parameter to control noise level. Default: seq(0.1, 1, 0.1)
#' @param type "NB", "N" or "LN". Noise type. "NB" means Negative binomial model,
#' "N" means normal model, "LN" means Log-normal model.
#' Default: "NB"
#'
#' @return  All the information will be written in the output path, and each file
#' is the generated data in a specific noise level.
#'
#' @importFrom tools file_path_sans_ext file_path_as_absolute
#' @importFrom progress progress_bar
#'
#' @export
#'
#' @examples
#' res <- pseudoExpr()
#' write.csv(x = res$mix, file = "mix.csv", row.names = T, quote = F)
#' addNoiseExpr(exprFile = "mix.csv", Pt = seq(0.1, 1, 0.1), type = "NB")
#'
addNoiseExpr <- function(exprFile,
                         outputPath = NULL,
                         prefix = NULL,
                         Pt = seq(0.1, 1, 0.1),
                         type = "NB") {
    if (is.null(prefix)) {
        prefix <- file_path_sans_ext(basename(exprFile))
        mess <- paste("The prefix will be", prefix, sep = " ")
        writeLines(mess)
    }

    if (is.null(outputPath)) {
        writeLines("Output path is not specified!")
        outputPath <- file.path(dirname(file_path_as_absolute(exprFile)),
                                prefix)
        mess <- paste("The output path will be:\n", outputPath, sep = "")
        writeLines(mess)
    }

    if (!file.exists(outputPath)) {
        dir.create(outputPath)
    }

    writeLines("Reading the expected expression value......")
    data <- read.csv(file = exprFile, header = TRUE, row.names = 1, encoding = "UTF-8")
    n_gene <- nrow(data)
    n_sample <- ncol(data)
    sample_names <- colnames(data)

    mess <- paste("There are", n_sample, "samples and", n_gene, "genes.", sep = " ")
    writeLines(mess)

    # output noise level 0 file (origin file)
    write.csv(x = data,
              file = file.path(outputPath, paste0(prefix, "_NL_0.csv")),
              row.names = TRUE)

    for (pt in Pt) {
        mess <- paste0("Generating noised data for pt = ", pt, "......")
        writeLines(mess)
        pb <- progress_bar$new(total = n_sample)
        this_data_noise <- list()
        for(i in seq(n_sample)) {
            pb$tick()
            this_data_noise[[sample_names[i]]] <- addNoise(x = data[[sample_names[i]]],
                                                           pt = pt,
                                                           type = type)
        }

        noisedData <- do.call(cbind, this_data_noise)
        rownames(noisedData) <- rownames(data)
        outputFileName <- paste0(prefix, "_NL_", pt, ".csv")
        write.csv(x = noisedData, file = file.path(outputPath, outputFileName), row.names = TRUE)
    }

}








