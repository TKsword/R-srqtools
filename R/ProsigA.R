#' @title Proteome Significance A
#'
#' When there is no biological duplication (or only 2 biological duplications), calculate the p-value of the proteomic data.
#'
#' Note: use empirical distribution.
#'
#' @param x A data.frame or matrix of proteomic data.
#' @param n Biological duplications number.
#'
#' @return P-value between each protein.
#' @export
#'
#' @examples

ProsigA <- function(x,n){
  df <- as.data.frame(x)
  if(n == 1){
    if(ncol(df) == 2){
      ratio <- df[,1] / df[,2]
    }else{
      return(warning("Check whether the number of input data is correct."))
    }
  }else if(n == 2){
    if(ncol(df) == 4){
      ratio <- (df[,1] + df[,2]) / (df[,3] + df[,4])
    }else{
      return(warning("Check whether the number of input data is correct."))
    }
  }else{
    return(warning("Enter the correct number of biological duplicates."))
  }
  ratio <- log2(ratio)
  order_ratio <- ratio[order(ratio)]
  quantiletmp <- quantile(order_ratio, c(0.1587,0.5,0.8413))
  rl <- as.numeric(quantiletmp[1])
  rm <- as.numeric(quantiletmp[2])
  rh <- as.numeric(quantiletmp[3])
  p <- unlist(lapply(ratio, function(x){
    if (x > rm){
      z <- (x-rm)/(rh-rm)
      pnorm(z,lower.tail = F)
    }else{
      z <- (rm-x)/(rm-rl)
      pnorm(z,lower.tail = F)
    }
  }))
  result <- data.frame(x,"P_value" = p)
  return(result)
}

