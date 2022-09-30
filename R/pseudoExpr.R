#' @title Generate pseudo mixture expression matrix
#'
#' @description Generate pseudo mixture expression matrix based on uniform
#' distribution, without any noise.
#'
#' @param n_sample Sample number to be generated, default: 50
#' @param n_gene Gene number to be generated, default: 1000
#' @param n_ct Cell type number to be generated, default: 10
#'
#' @return a list, prop means the pseudo proportion, ref means the pseudo
#' external reference, mix is the output of ref %*% prop
#'
#' @export
#'
#' @examples
#' res <- pseudoExpr()
#'
pseudoExpr <- function(n_sample = 50, n_gene = 1000, n_ct = 10) {

    ref <- matrix(data = sample(x = 10000, size = n_gene * n_ct, replace = TRUE),
                  nrow = n_gene,
                  ncol = n_ct)
    rownames(ref) <- paste0("gene", seq(nrow(ref)))
    colnames(ref) <- paste0("ct", seq(ncol(ref)))

    prop <- matrix(data = sample(x = 1000, size = n_sample * n_ct, replace = TRUE),
                   nrow = n_ct, ncol = n_sample)
    prop <- apply(X = prop, MARGIN = 2, FUN = v_norm)
    rownames(prop) <- paste0("ct", seq(nrow(prop)))
    colnames(prop) <- paste0("S", seq(ncol(prop)))

    mix <- ref %*% prop

    return(list(mix = mix, ref = ref, prop = prop))
}











