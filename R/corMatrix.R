#' @title Calculate the correlation coefficient between two Matrices
#'
#' @description To calculate the correlation coefficient between two equivalent matrices or non-equivalent matrices, refer to the specific method:
#' @description Indahl, U.G.; Næs, T.; Liland, K.H.; 2018. A similarity index for comparing coupled matrices. Journal of Chemometrics; e3049.
#'
#' @param X1 matrix 1 (data.frame format is also accepted)
#' @param X2 matrix 2 (data.frame format is also accepted)
#' @param equal whether the two matrices are consistent in direction, that is, the equivalent matrix.
#' @param method method of calculating correlation coefficient: "SMI", "RV", "RV2", "RVadj", "PSI", "pearson", "r1", "r2", "r3", "r4", "spearman", "GCD" or "kendall"
#'
#' @param ncomp1 maximum number of subspace components from the first matrix.
#' @param ncomp1 maximum number of subspace components from the second matrix.
#'
#'
#'
#' @export corMatrix
#'
#' @import dplyr
#' @import stats
#' @import MatrixCorrelation
#'
#'
#' @examples

corMatrix <- function(X1, X2, equal = TRUE, method = "pearson",
                      ncomp1, ncomp2, ...){
  if (equal == FALSE){
    if (missing(ncomp1) && ("SMI" %in% methods || "GCD" %in% methods)) {
      ncomp1 <- min(5, Rank(X1) - 1)
      warning(paste("ncomp1 not specified, defaulting to ", ncomp1, sep = ""))
    }
    if (missing(ncomp2) && ("SMI" %in% methods || "GCD" %in% methods)) {
      ncomp2 <- min(5, Rank(X2) - 1)
      warning(paste("ncomp2 not specified, defaulting to ", ncomp2, sep = ""))
    }
    n.met <- length(methods)
    results <- numeric(n.met)
    for (i in 1:n.met) {
      results[i] <- switch(methods[i], SMI = SMI(X1, X2, ncomp1, ncomp2, ...)[ncomp1, ncomp2],
                           RV = RV(X1, X2),
                           RV2 = RV2(X1, X2),
                           RVadj = RVadj(X1, X2),
                           PSI = PSI(X1, X2),
                           r1 = r1(X1, X2),
                           r2 = r2(X1, X2),
                           r3 = r3(X1, X2),
                           r4 = r4(X1, X2),
                           GCD = GCD(X1, X2, ncomp1, ncomp2))
    }
    names(results) <- methods
    if (!is.null(digits) && !is.na(digits)) {
      results <- round(results, digits)
    }
  }else{
    X1 <- as.matrix(X1) %>% as.vector(.)
    X2 <- as.matrix(X2) %>% as.vector(.)
    if (length(X1) != length(X2)){
      return(errorCondition("Please check if the two matrices are equivalent"))
    }
    result <- cor(X1, X2, method = method)
  }
  return(result)
}
