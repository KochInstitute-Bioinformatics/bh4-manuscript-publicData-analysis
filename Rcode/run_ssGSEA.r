
ssGSEA.project.dataset(input.ds = "XXX.gct", output.ds = "customSets_ssGSEA", gene.sets.dbfile.list = "customSets.gmx")

rc.toPlot <- calculate_zscore(toPlot)
round(mean(rowMeans(rc.toPlot)),2)
round(mean(rowVars(rc.toPlot)),2)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")