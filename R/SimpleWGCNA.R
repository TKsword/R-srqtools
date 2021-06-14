#' A quick implementation function of WGCNA
#'
#' The purpose of this function is to calculate the result of WGCNA simply and quickly based on the expression data. Most of the parameters refer to the original document of WGCNA, and a few parameters are set according to the author's own experience.
#'
#' Note: This function cannot replace the complete WGCNA solution. If the purpose is not to get a quick overview but to draw a formal conclusion, please consider its applicability to WGCNA.
#'
#' @param dataExpr A data.frame of gene expression matrix
#' @param PNcor Logical value, whether the positive and negative correlations are considered separately
#' @param useEmpiricalBeta Logical value, whether to use empirical beta value
#' @param FileName The name of the output file
#'
#' @return
#' @export
#'
#' @import WGCNA
#' @import stringr
#' @import openxlsx
#'
#'
#' @examples
SimpleWGCNA <- function(dataExpr, PNcor = TRUE, useEmpiricalBeta = FALSE, FileName = "SimpleWGCNA result"){
  nSamples <- ncol(dataExpr)
  if(nSamples < 10){
    return(warning("Please use more samples to increase diversity and improve the accuracy of analysis."))
  }
  type <- ifelse(PNcor == TRUE, "signed", "unsigned")
  options(stringsAsFactors = FALSE)
  dataExpr <- as.data.frame(t(dataExpr))
  gsg <- goodSamplesGenes(dataExpr, verbose = 3)
  if(!gsg$allOK){
    if (sum(!gsg$goodGenes) > 0)
      printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples) > 0)
      printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", ")));
    dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  sampleTree <- hclust(dist(dataExpr), method = "average")
  cutHeight <- max(sampleTree$height) * 0.7
  clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 2)
  keepSamples <- (clust >= 1)
  dataExpr <- dataExpr[keepSamples, ]
  powers <- c(1:20)
  sft <- pickSoftThreshold(dataExpr, powerVector = powers, networkType = type, verbose = 5)
  fitIndices <- sft$fitIndices
  if(useEmpiricalBeta == FALSE){
    meank <- fitIndices[fitIndices$mean.k. >= 100,]
    power <- meank[(-sign(meank[,3])*meank[,2]) == abs(max(-sign(meank[,3])*meank[,2])),1]
    power <- ifelse(length(power) == 0, 0, power)
    if(power == 0){
      return(warning("Unable to calculate the beta value that conforms to the non-scale network distribution, please reorganize the data, or try to use the empirical beta value."))
    }else if(fitIndices$SFT.R.sq[power] < 0.5){
      return(warning("Unable to calculate the beta value that conforms to the non-scale network distribution, please reorganize the data, or try to use the empirical beta value."))
    }else{}
  }else{
    power <- ifelse(nSamples < 20, ifelse(type == "unsigned", 9, 18),ifelse(nSamples < 30, ifelse(type == "unsigned", 8, 16),
                                                                            ifelse(nSamples < 40, ifelse(type == "unsigned", 7, 14),
                                                                                   ifelse(type == "unsigned", 6, 12))))
  }
  net <- blockwiseModules(dataExpr, power = power,
                          TOMType = type, minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, corType = "pearson",
                          verbose = 3)
  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  pdf(file = "Cluster Dendrogram.pdf")
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  MEs <- net$MEs
  MEs_col <- MEs
  colnames(MEs_col) <- paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  if(ncol(MEs_col) == 3){
    pdf(file = "Eigengene adjacency heatmap.pdf")
    plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                          marDendro = c(3,3,2,4),
                          marHeatmap = c(3,4,2,2), plotDendrograms = F,
                          xLabelsAngle = 90)
    dev.off()
  }else if(ncol(MEs_col) > 3){
    pdf(file = "Eigengene adjacency heatmap.pdf")
    plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                          marDendro = c(3,3,2,4),
                          marHeatmap = c(3,4,2,2), plotDendrograms = T,
                          xLabelsAngle = 90)
    dev.off()
  }else{
    warning("Insufficient number of classified modules.")
  }
  KME <- signedKME(dataExpr, MEs_col, outputColumnName = "")
  colorNames <- table(moduleColors)
  result <- list()
  for(i in 1:length(names(colorNames))){
    cME <- KME[moduleColors == names(colorNames)[i],colnames(KME) == names(colorNames)[i]]
    cME <- cbind(rownames(KME[moduleColors == names(colorNames)[i],]), cME)
    colnames(cME) <- c("GeneSymbol", names(colorNames)[i])
    result[[i]] <- cME
  }
  write.xlsx(result,paste0(FileName,".xlsx"))
  return()
}
