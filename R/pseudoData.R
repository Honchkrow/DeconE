#' @title Generate pseudo data for plotting
#'
#' @description Generate pseudo data (prediction values) for plotting
#'
#' @param type 1, 2 or 3 for different usage
#' @param outputPath path to save tmp file
#'
#' @return The pseudo data
#'
#' @export
#'
#' @examples
#' res <- pseudoData(type = 1)
#'
pseudoData <- function(type = 1, outputPath = 'test_file') {
    if (!dir.exists(outputPath)) {
        dir.create(outputPath)
    }

    if (type == 1) {
        actual <- matrix(data = seq(50), nrow = 5)
        colnames(actual) <- paste0("S", seq(10))
        rownames(actual) <- paste0("ct", seq(5))

        predicted <- matrix(data = seq(50), nrow = 5) + 1
        colnames(predicted) <- paste0("S", seq(10))
        rownames(predicted) <- paste0("ct", seq(5))

        res <- list(actual = actual, predicted = predicted)
    } else if (type == 2) {
        actual <- matrix(data = sample(x = seq(1000), size = 50), nrow = 5)
        colnames(actual) <- paste0("S", seq(10))
        rownames(actual) <- paste0("ct", seq(5))

        write.csv(x = as.data.frame(actual),
                  file = file.path(outputPath, "actual.csv"),
                  row.names = TRUE)
        actual <- file.path(outputPath, "actual.csv")

        predicted <- c()
        for (i in seq(10)) {
            tmp_pred <- matrix(data = sample(x = seq(1000), size = 50), nrow = 5)
            colnames(tmp_pred) <- paste0("S", seq(10))
            rownames(tmp_pred) <- paste0("ct", seq(5))
            tmp_file <- file.path(outputPath, paste0("pred_", i, ".csv"))
            write.csv(x = as.data.frame(tmp_pred),
                      file = tmp_file,
                      row.names = TRUE)
            predicted <- c(predicted, tmp_file)
        }

        noise_level <- paste0("NL", seq(0.1, 1, 0.1))

        res <- list(actual = actual, predicted = predicted, noise_level = noise_level)
    } else if (type == 3) {
        actual <- matrix(data = sample(x = seq(1000), size = 50), nrow = 5)
        colnames(actual) <- paste0("S", seq(10))
        rownames(actual) <- paste0("ct", seq(5))

        write.csv(x = as.data.frame(actual),
                  file = file.path(outputPath, "actual.csv"),
                  row.names = TRUE)
        actual <- file.path(outputPath, "actual.csv")

        predicted <- list()
        for (m in c("method1", "method2", "method3", "method4", "method5", "method6")) {
            tmp_v <- c()
            for (i in seq(10)) {
                tmp_pred <- matrix(data = sample(x = seq(1000), size = 50), nrow = 5)
                colnames(tmp_pred) <- paste0("S", seq(10))
                rownames(tmp_pred) <- paste0("ct", seq(5))
                tmp_file <- file.path(outputPath, paste0(m, "_", i, ".csv"))
                write.csv(x = as.data.frame(tmp_pred),
                          file = tmp_file,
                          row.names = TRUE)
                tmp_v <- c(tmp_v, tmp_file)
            }
            predicted[[m]] <- tmp_v
        }
        noise_level <- paste("NL", seq(1, 10), sep = "")
        res <- list(actual = actual, predicted = predicted, noise_level = noise_level)
    }

    return(res)
}
