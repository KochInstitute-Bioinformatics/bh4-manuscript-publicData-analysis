KSScore = function(expMat, gseID, outDirectory) {

  cat("Calculating EMT Score by KS score method ....\n")
  
  # Remove rows with NaN or -Inf
  remIdx = which(apply(expMat, 1, function(x) any(is.nan(x) | x == -Inf)))
  if (length(remIdx) > 0) expMat = expMat[-remIdx, ]
  
  sampleNum = ncol(expMat)
  genes = rownames(expMat)
  exp = apply(expMat[, 1:sampleNum], 2, as.numeric)
  
  EMTSignature = data.frame(read.xlsx(KSScore.EMTSignature.file, colNames = FALSE))
  commonSig = intersect(EMTSignature[, 1], genes)
  EMTExpIdx = match(commonSig, genes)
  EMTExp = exp[EMTExpIdx, ]
  print(dim(EMTExp))
  EMTGrpIdx = match(commonSig, EMTSignature[, 1])
  geneCat = EMTSignature[EMTGrpIdx, 2]
  epiIdx = which(geneCat == "Epi")
  mesIdx = which(geneCat == "Mes")
  
  message(length(commonSig), " KS score genes in our expression matrix.")
  
  # Perform KS test
  sampleScore2 = matrix(0, nrow = ncol(EMTExp), ncol = 6)
  rownames(sampleScore2) = colnames(EMTExp)
  
  for (i in 1:ncol(EMTExp)) {
    ksTwoSided = ks.test(EMTExp[mesIdx, i], EMTExp[epiIdx, i])
    ksResGrt = ks.test(EMTExp[mesIdx, i], EMTExp[epiIdx, i], alternative = "greater")
    ksResLess = ks.test(EMTExp[epiIdx, i], EMTExp[mesIdx, i], alternative = "greater")
    
    sampleScore2[i, ] = c(
      ksTwoSided$statistic, ksTwoSided$p.value,
      ksResGrt$statistic, ksResGrt$p.value,
      ksResLess$statistic, ksResLess$p.value
    )
  }
  
  # Assign signs to EMT score
  finalScore = numeric(nrow(sampleScore2))
  names(finalScore) = rownames(sampleScore2)
  
  for (i in 1:nrow(sampleScore2)) {
    if (sampleScore2[i, 4] < 0.05) {
      finalScore[i] = -1 * sampleScore2[i, 3]
    } else if (sampleScore2[i, 6] < 0.05) {
      finalScore[i] = sampleScore2[i, 5]
    } else {
      maxVal = max(sampleScore2[i, 3], sampleScore2[i, 5])
      finalScore[i] = ifelse(sampleScore2[i, 5] == maxVal, maxVal, -1 * maxVal)
    }
  }
  
  # Write to file
  outputFile = file.path(outDirectory, paste0(gseID, "_EMT_KS.txt"))
  write.table(
    cbind(Sample = names(finalScore), EMT_KS_Score = finalScore),
    file = outputFile, sep = '\t', row.names = FALSE, quote = FALSE
  )
  
  # Return named vector of scores
  return(finalScore)
}
