rm(list = ls())
load("TCGA_OV_count_transformation.Rdata")
write.table(vsdd, "TCGA_OV_count_transformation.txt")
