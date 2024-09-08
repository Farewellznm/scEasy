dir.create("examples")

library(scEasy)



gl = read.csv("../../datebase/genesets/aging_genesets/Aging Gene Sets.csv")
gl = gl$Symbol

scEasy::Run_Gene_Enrichment(geneList = gl[1:50],org = "hsa")
