#' @title Add noise to gene expression data
#'
#' @description Add noise based on negtive bionomial distribution. Please see
#' "A benchmark for RNA-seq deconvolution analysis
#' under dynamic testing environments" in Genome Biology.
#'
#' @param x a gene expression numeric vector.
#' @param pt parameter to control noised level. Default: 0.1
#' @param type "NB", "N" or "LN". "NB" means Negative binomial model,
#' "N" means normal model, "LN" means Log-normal model.
#' Default: "NB"
#'
#' @return the sample length vector.
#'
#' @importFrom stats rnorm rpois
#'
#' @export
#'
#' @examples
#' res <- addNoise(x = seq(100))
#'
addNoise <- function (x = NULL, pt = 0.1, type = "NB") {
    x_size <- length(x)

    if (type == "NB") {
        # step1 compute sigma
        sigma_tmp <- 1.8 * pt + (1 / sqrt(x))
        delta <- rnorm(n = x_size, mean = 0, sd = 0.25)
        sigma <- sigma_tmp * exp(delta / 2)

        # step2 compute u
        shape <- 1 / (sigma ** 2)
        scale <- x / shape
        u <- mapply(FUN = G, shape = shape, scale = scale)
        u[is.na(u)] <- 0

        # step3 generate v
        v <- sapply(X = u, FUN = rpois, n = 1)
    } else if (type == "N") {
        # step1 compute normal
        n1 <- rnorm(n = x_size, mean = 0, sd = sqrt(10 * pt))

        # step2
        n2 <- log2((x + 1)) + n1
        v <- 2 ** n2
    } else if (type == "LN") {
        # step1
        n1 <- rnorm(n = x_size, mean = 0, sd = sqrt(10 * pt))
        v <- x + 2 ** n1

    } else {
        stop("Parameter 'type' is invalid!")
    }


    return(v)
}


#' @title Generate a gamma random number
#'
#' @description Generate a gamma random number.
#'
#' @param shape shape
#' @param scale scale
#'
#' @return random number from gamma distribution
#'
#' @importFrom stats rgamma
#'
#' @keywords internal
#'
G <- function (shape, scale) {
    out <- rgamma(n = 1, shape = shape, scale = scale)
    return(out)
}


#' @title Normalize a vector
#'
#' @description Normalize a vector
#'
#' @param x a numeric vector
#'
#' @return Normalized numeric vector
#'
#' @export
#'
#' @examples
#' res <- v_norm(seq(100))
#'
v_norm <- function(x = NULL){
    normRes <- x / sum(x)
    return(normRes)
}


#' @title Merge dataframes
#'
#' @description Merge and align dataframes.
#'
#' @param x A dataframe
#' @param ... Other dataframes
#' @param by Flag used to merge, default: "row.names"
#'
#' @return A dataframe
#'
#' @keywords internal
#'
merge.all <- function(x, ..., by = "row.names"){
    L <- list(...)
    for (i in seq_along(L)) {
        x <- merge(x, L[[i]], by = by)
        rownames(x) <- x$Row.names
        x$Row.names <- NULL
    }
    return(x)
}


#' @title Convert a gene expression count matrix into TPM matrix.
#'
#' @description Convert a gene expression count matrix into TPM matrix
#'
#' @param data A dataframe or matrix object.
#' Note: Must contain a column named "Length", "Length" means gene length.
#'
#' @return A dataframe or matrix object.
#'
#' @export
#'
#' @examples
#' data <- data.frame(s1 = seq(100), Length = seq(100))
#' res <- TPM(data)
#'
TPM <- function(data = NULL){
    if("Length" %in% colnames(data)){
        rawCounts <- within(data, rm("Length"))
        length <- data["Length"]
        normCounts <- rawCounts / c(length)
        normCounts <- apply(X = normCounts, MARGIN = 2, FUN = v_norm) * 1e6
    } else {
        stop("Must contain a column named 'Length', 'Length' means gene length.")
    }

    return(normCounts)
}


#' @title Convert a gene expression count matrix into CPM matrix.
#'
#' @description Convert a gene expression count matrix into CPM matrix
#'
#' @param data A dataframe or matrix object.
#' Note: Must contain a column named "Length", "Length" means gene length.
#'
#' @return A dataframe or matrix object.
#'
#' @export
#'
#' @examples
#' data <- data.frame(s1 = seq(100), Length = seq(100))
#' res <- CPM(data)
#'
CPM <- function(data = NULL){
    if("Length" %in% colnames(data)){
        rawCounts <- within(data, rm("Length"))
    }else{
        stop("Must contain a column named 'Length', 'Length' means gene length.")
    }
    normCounts <- apply(X = rawCounts, MARGIN = 2, FUN = v_norm) * 1e6
    return(normCounts)
}


#' @title Split a dataframe into train and test set.
#'
#' @description Split a dataframe into train and test set based on the labelCol,
#' 1 sample in one class at least.
#'
#' @param df A dataframe for all cell type data.
#' @param p Proportion of sample in train set.
#' @param labelCol column name to record cell types.
#' @param outVarCol column name to output, usually is sample name like RUN Accession Number.
#'
#' @return A list contains train and sample names.
#'
#' @keywords internal
#'
multiClassSample <- function(df = NULL, p = NULL,
                             labelCol = "Required_Cell_Type",
                             outVarCol = "RUN_Accession"){
    df.summary <- data.frame(table(df[[labelCol]]))
    df.summary[["train"]] <- as.integer(df.summary[["Freq"]] * p)
    df.summary[["test"]] <- df.summary[["Freq"]] - df.summary[["train"]]
    classNames <- unique(df[[labelCol]])

    train <- c()
    test <- c()

    for(i in seq(nrow(df.summary))){
        samples <- df[df[[labelCol]] == df.summary$Var1[i], ]
        train.samples <- sample(x = samples[[outVarCol]], size = df.summary$train[i])
        train <- c(train, train.samples)
        test <- c(test, setdiff(samples[[outVarCol]], train.samples))
    }

    return(list(train = train, test = test))
}



#' @title Generate reference matrix.
#'
#' @description Generate reference matrix for deconvolution. This function computes
#' the expression and expression variance matrix for different cell types.
#'
#' @param df A dataframe for cell type expression data.
#' @param attr Attribute information.
#' @param cellTypes Cell types information.
#'
#' @return A reference list contains expression and expression variance matrix.
#'
#' @importFrom stats sd
#'
#' @keywords internal
#'
refGenerator <- function(df = NULL, attr = NULL, cellTypes = NULL){
    ref <- as.data.frame(matrix(nrow = nrow(df), ncol = 0))
    rownames(ref) <- rownames(df)

    std <- as.data.frame(matrix(nrow = nrow(df), ncol = 0))
    rownames(std) <- rownames(df)

    for(this.celltype in cellTypes){
        this.SRR <- attr$RUN_Accession[which(attr$Required_Cell_Type == this.celltype)]
        this.selected <- this.SRR[which(this.SRR %in% colnames(df))]
        ref[this.celltype] <- apply(X = df[, this.selected], MARGIN = 1, FUN = mean)
        std[this.celltype] <- apply(X = df[, this.selected], MARGIN = 1, FUN = sd)
    }

    ref <- as.matrix(ref)
    std <- as.matrix(std)

    return(list(ref = ref, std = std))
}









