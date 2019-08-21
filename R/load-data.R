# where is everything?
setwd("~/Dropbox/research/functional-gwas/")

# load packages
library(lubridate)

# load raw data sets
age_data <- read.csv("data/raw/GP_2018_covariates_wAgeData_updated_wsibgroups_qvalues_inddiv_approx_coords191018.csv",
                     stringsAsFactors = FALSE,
                     row.names = 1)
snp_data <- read.table("data/raw/gp_snp_data_012_format_3775loci.txt",
                       stringsAsFactors = FALSE,
                       row.names = 1,
                       sep = " ")

# match rows of data sets
growth_data <- age_data[, grep("ZL|Primordium", colnames(age_data))]
growth_measured <- !is.na(growth_data$Primordium)
growth_data <- growth_data[growth_measured, ]
rownames(growth_data) <- age_data$id[growth_measured]
colnames(growth_data) <- c(paste("Y", seq_len(ncol(growth_data) - 1), sep = "_"))
full_columns <- apply(growth_data, 2, function(x) sum(is.na(x))) < nrow(growth_data)
growth_data <- growth_data[, full_columns]

# create a covariates file
covar_names <- c("id", "pop", "sex", "oto_age",
                 "site", "latitude", "longitude",
                 "date_start",
                 "total_length_mm", "weight_g")
covariate_data <- age_data[, covar_names]
covariate_data <- covariate_data[growth_measured, ]
covariate_data$date_formatted <- parse_date_time(covariate_data$date_start, orders = c("dmy"))
covariate_data$date_start <- NULL

# pull out summary genetics info
genetic_data <- age_data[, grep("^Q|PHt|Hs|IR|HL", colnames(age_data))]
genetic_data <- genetic_data[growth_measured, ]

# filter down to those individuals with growth data as well
snp_data <- snp_data[match(rownames(growth_data), rownames(snp_data)), ]
to_keep <- !is.na(snp_data[, 1])
growth_data <- growth_data[to_keep, ]
covariate_data <- covariate_data[to_keep, ]
genetic_data <- genetic_data[to_keep, ]
snp_data <- snp_data[to_keep, ]

# write output files
write.csv(growth_data, file = "data/compiled/gp_fgwas_pheno.csv", row.names = FALSE)
write.csv(snp_data, file = "data/compiled/gp_fgwas_snps.csv", row.names = FALSE)
write.csv(covariate_data, file = "data/compiled/gp_fgwas_covariates.csv", row.names = FALSE)
write.csv(genetic_data, file = "data/compiled/gp_fgwas_genetics.csv", row.names = FALSE)
