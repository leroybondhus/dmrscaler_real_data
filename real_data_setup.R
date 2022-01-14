## real_data_setup.R ##
# perform preprocessing on all datasets and set up method sets


## install DMRscaler
library("devtools")
library("roxygen2")
install("../DMRscaler", quick=T)

library(doParallel)
registerDoParallel()

results_dir <-paste("./results/")


## Example Dataset Setup
##We will use data from GSE149960 from (REFERENCE) with DNA methylation from fibroblasts from progeria patients and controls measured on the Illumina methylation EPIC array.

library("GEOquery")



# Pre-processing
###  Reading of idat files done with minfi library ###

gse <- getGEO("GSE74432", GSEMatrix = TRUE)
phen <- gse$GSE74432_series_matrix.txt.gz@phenoData@data
phen <- phen[intersect( grep("[Cc]ontrol", phen$characteristics_ch1.1),grep("whole blood",phen$characteristics_ch1.2) ),]
rm(gse)

## get methylation data as idat files (NOTE: this saves files locally in working directory, unpacked size is 2.01 Gb
if(length(list.files("GSE74432/idat", pattern = "idat$"))==0){
  getGEOSuppFiles("GSE74432")
  untar("GSE74432/GSE74432_RAW.tar", exdir = "GSE74432/idat")
  file.remove("GSE74432/GSE74432_RAW.tar")
  idat_files <- list.files("GSE74432/idat", pattern = "idat.gz$", full = TRUE)
  sapply(idat_files, gunzip, overwrite = TRUE); rm(idat_files)
}

##  Preprocessing
library("minfi")


###  Reading of idat files done with minfi library ###
idats_dir<-"GSE74432/idat"
targets <- data.frame("Basename"= stringr::str_split_fixed(basename(phen$supplementary_file), "_Grn", 2)[,1] )
RGSet <- read.metharray.exp(base = idats_dir, targets = targets)
GRset.funnorm <- preprocessFunnorm(RGSet);rm(RGSet)
snps <- getSnpInfo(object = GRset.funnorm)
GRset.funnorm <- dropLociWithSnps(GRset.funnorm, snps=c("SBE", "CpG"), maf=0);rm(snps)
rm(idats_dir)


controls <- grep("[Cc]ontrol",phen$title)
locs <- getLocations(GRset.funnorm)
locs <- data.frame("names"=locs@ranges@NAMES, "pos"=locs@ranges@start, "chr" = rep(locs@seqnames@values, locs@seqnames@lengths))
B <- getBeta(GRset.funnorm)


method_set_list <- list(
  dmrscaler_1 = list(method="dmrscaler",
                     function_call=" DMRscaler::dmrscaler(locs = locs, locs_pval_cutoff = 0.01, region_signif_method = \"bon\", region_signif_cutoff = 0.01, window_type = \"k_near\", window_sizes = c(2,4,8,16,32,64), output_type = \"comp\") " ),
  dmrscaler_2 = list(method="dmrscaler",
                     function_call=" DMRscaler::dmrscaler(locs = locs, locs_pval_cutoff = 0.01, region_signif_method = \"ben\", region_signif_cutoff = 0.01, window_type = \"k_near\", window_sizes = c(2,4,8,16,32,64), output_type = \"comp\") " ),
  bumphunter_1 = list(method="bumphunter",
                      function_call="bumphunter(B_mod,as.matrix(design),chr = locs$chr, pos=locs$pos, cutoff=0.1, maxGap=1e3, B=250, smoothFunction=loessByCluster )"),
  bumphunter_2 = list(method="bumphunter",
                      function_call="bumphunter(B_mod,as.matrix(design),chr = locs$chr, pos=locs$pos, cutoff=0.1, maxGap=1e6, B=250, smoothFunction=loessByCluster )"),
  dmrcate_1 = list(method="dmrcate",
                   function_call="dmrcate(myannotation, lambda=1e3, C=2)"),
  dmrcate_2 = list(method="dmrcate",
                   function_call="dmrcate(myannotation, lambda=1e6, C=2000)"),
  combp_1 = list(method="combp",
                 function_call="system(\"comb-p pipeline -c 4 --dist 1000 --step 100 --seed 1e-3 --region-filter-p 0.1 -p " ),
  combp_2 = list(method="combp",
                 function_call="system(\"comb-p pipeline -c 4 --dist 10000 --step 5000 --seed 1e-3 --region-filter-p 0.1 -p " )
)


## write Rdata objects to load into each real_data_individual_run.R script
filename <- paste("real_data_setup.Rdata")
save.image(file=filename)
## write dataset_table
dataset_table <- names(dataset_list)
filename <- paste("dataset_table.csv")
write.table(dataset_table, filename, row.names = F, col.names=F, sep=",")

## write method_table
method_table <- names(method_set_list)
filename <- paste("method_table.csv")
write.table(method_table, filename, row.names = F, col.names=F, sep=",")
