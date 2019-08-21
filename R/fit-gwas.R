# where am I working?
setwd("~/Dropbox/research/functional-gwas/")

# load the fgwas package indirectly (many R scripts)
file_list <- dir("~/Downloads/fgwas/R/")
for (i in seq_along(file_list))
  source(paste0("~/Downloads/fgwas/R/", file_list[i]))

# load some other packages we need
library(mvtnorm)
library(MSBVAR)

# run GLS model
GLSApp.onestep("data/compiled/gp_fgwas.tped",
               "data/compiled/gp_fgwas.tfam",
               "data/compiled/gp_fgwas_pheno.csv",
               file.rdata = "outputs/test-",
               fields = list(Y = "ZL"))
