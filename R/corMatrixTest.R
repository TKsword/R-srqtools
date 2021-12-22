#' @title Hypothesis test on the correlation number between two matrices
#'
#'
#' @param X1 matrix 1 (data.frame format is also accepted)
#' @param X2 matrix 2 (data.frame format is also accepted)
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. "greater" corresponds to positive association, "less" to negative association.
#' @param conf.level confidence level for the returned confidence interval. Currently only used for the Pearson product moment correlation coefficient if there are at least 4 complete pairs of observations.
#'
#'
#'
#' @export corMatrixTest
#'
#' @examples

corMatrixTest <- function(X1, X2,
                          alternative = "two.sided",
                          method = "pearson",
                          exact = NULL,
                          conf.level = 0.95, ...){
  X1 <- as.matrix(X1) %>% as.vector(.)
  X2 <- as.matrix(X2) %>% as.vector(.)
  if (length(X1) != length(X2)){
    return(errorCondition("Please check if the two matrices are equivalent"))
  }else{
    result.full <- cor.test(X1, X2, alternative = alternative, method = method, exact = exact, conf.level = conf.level, ...)
  }
  result <- c(round(result.full$p.value, 5), round(result.full$conf.int, 5))
  names(result) <- c("P-value", "conf.int.lower", "conf.int.higher")
  return(result)
}
