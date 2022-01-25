## real_data_setup.R ##
# perform preprocessing on all datasets and set up method sets


## install DMRscaler
library("devtools")
library("roxygen2")
install("../DMRscaler", quick=T)

library(doParallel)
registerDoParallel()

results_dir <-paste("./results/")

data_set_list <- list()

library("GEOquery")


# Pre-processing
###  Reading of idat files done with minfi library ###

gse <- getGEO("GSE74432", GSEMatrix = TRUE)
phen <- gse$GSE74432_series_matrix.txt.gz@phenoData@data
phen <- phen[grep("whole blood",phen$characteristics_ch1.2),]
phen$col_names <- stringr::str_split_fixed(basename(phen$supplementary_file), "_Grn", 2)[,1]
rm(gse)

## get methylation data as idat files (NOTE: this saves files locally in working directory, unpacked size is 2.01 Gb
if(length(list.files("../dmrscaler_simulation/GSE74432/idat", pattern = "idat$"))==0){
  getGEOSuppFiles("GSE74432")
  untar("GSE74432/GSE74432_RAW.tar", exdir = "GSE74432/idat")
  file.remove("GSE74432/GSE74432_RAW.tar")
  idat_files <- list.files("GSE74432/idat", pattern = "idat.gz$", full = TRUE)
  sapply(idat_files, gunzip, overwrite = TRUE); rm(idat_files)
}

##  Preprocessing
library("minfi")

###  Reading of idat files done with minfi library ###
idats_dir<-"../dmrscaler_simulation/GSE74432/idat"

####   Sotos
phen_subset <- phen[union(grep("Sotos", phen$`disease state:ch1`),
                            grep("Control", phen$`disease state:ch1`)),]
targets <- data.frame("Basename"= phen_subset$col_names)
RGSet <- read.metharray.exp(base = idats_dir, targets = targets)
GRset.funnorm <- preprocessFunnorm(RGSet);rm(RGSet)
snps <- getSnpInfo(object = GRset.funnorm)
GRset.funnorm <- dropLociWithSnps(GRset.funnorm, snps=c("SBE", "CpG"), maf=0);rm(snps)
g1 <- phen_subset$col_names[grep("[Cc]ontrol", phen_subset$`disease state:ch1`)]
g2 <- phen_subset$col_names[grep("[Ss]otos", phen_subset$`disease state:ch1`)]
locs <- getLocations(GRset.funnorm)
locs <- data.frame("names"=locs@ranges@NAMES, "pos"=locs@ranges@start, "chr" = rep(locs@seqnames@values, locs@seqnames@lengths))
B <- getBeta(GRset.funnorm)

data_set_list[["Sotos"]] <- list( data_set_name = "Sotos", g1 = g1, g2 = g2, locs = locs, B = B,
                                  g1g2_labels = data.frame(g1="Control",g2="Sotos"))

####   Weaver
phen_subset <- phen[union(grep("Weaver", phen$`disease state:ch1`),
                          grep("Control", phen$`disease state:ch1`)),]
targets <- data.frame("Basename"= phen_subset$col_names)
RGSet <- read.metharray.exp(base = idats_dir, targets = targets)
GRset.funnorm <- preprocessFunnorm(RGSet);rm(RGSet)
snps <- getSnpInfo(object = GRset.funnorm)
GRset.funnorm <- dropLociWithSnps(GRset.funnorm, snps=c("SBE", "CpG"), maf=0);rm(snps)
g1 <- phen_subset$col_names[grep("[Cc]ontrol", phen_subset$`disease state:ch1`)]
g2 <- phen_subset$col_names[grep("[Ww]eaver", phen_subset$`disease state:ch1`)]
locs <- getLocations(GRset.funnorm)
locs <- data.frame("names"=locs@ranges@NAMES, "pos"=locs@ranges@start, "chr" = rep(locs@seqnames@values, locs@seqnames@lengths))
B <- getBeta(GRset.funnorm)

data_set_list[["Weaver"]] <- list( data_set_name = "Weaver", g1 = g1, g2 = g2, locs = locs, B = B,
                                  g1g2_labels = data.frame(g1="Control",g2="Weaver"))

####   Sex
phen_subset <- phen[ grep("Control", phen$`disease state:ch1`),]
targets <- data.frame("Basename"= phen_subset$col_names)
RGSet <- read.metharray.exp(base = idats_dir, targets = targets)
GRset.funnorm <- preprocessFunnorm(RGSet);rm(RGSet)
snps <- getSnpInfo(object = GRset.funnorm)
GRset.funnorm <- dropLociWithSnps(GRset.funnorm, snps=c("SBE", "CpG"), maf=0);rm(snps)
g1 <- phen_subset$col_names[grep("[Ff]emale", phen_subset$`gender:ch1`)]
g2 <- phen_subset$col_names[grep("[Ff]emale", invert = T, phen_subset$`gender:ch1`)]
locs <- getLocations(GRset.funnorm)
locs <- data.frame("names"=locs@ranges@NAMES, "pos"=locs@ranges@start, "chr" = rep(locs@seqnames@values, locs@seqnames@lengths))
B <- getBeta(GRset.funnorm)

data_set_list[["Sex"]] <- list( data_set_name = "Sex", g1 = g1, g2 = g2, locs = locs, B = B,
                                   g1g2_labels = data.frame(g1="Female",g2="Male"))



####   KAT6A
idats_dir<-"/u/project/arboleda/DATA/DNAme_Array/KAT6A_DNA_methylation_data/idat_KAT6A_and_control"
idats_files<-list.files(path=idats_dir, pattern = "*.idat")
phen <- read.metharray.sheet("/u/project/arboleda/DATA/DNAme_Array/KAT6A_DNA_methylation_data/", pattern = "Formatted", verbose = TRUE)
phen$Basename<-paste(phen$Slide, "_", phen$Array, sep = "")
#Extract only KAT6A patients and Controls for analysis
phen <- phen[c(grep("^CONTROL",phen$Sample_Name),grep("^KAT6",phen$Sample_Name)),]
phen <- phen[-grep("-2", phen$Sample_Name),] ## remove technical replicates

targets <- data.frame("Basename"= phen$Basename)
RGSet <- read.metharray.exp(base = idats_dir, targets = targets)
GRset.funnorm <- preprocessFunnorm(RGSet);rm(RGSet)
snps <- getSnpInfo(object = GRset.funnorm)
GRset.funnorm <- dropLociWithSnps(GRset.funnorm, snps=c("SBE", "CpG"), maf=0);rm(snps)
g1 <- phen$Basename[grep("Control", phen$Sample_Name, ignore.case = T)]
g2 <- phen$Basename[grep("KAT6A", phen$Sample_Name, ignore.case = T)]
locs <- getLocations(GRset.funnorm)
locs <- data.frame("names"=locs@ranges@NAMES, "pos"=locs@ranges@start, "chr" = rep(locs@seqnames@values, locs@seqnames@lengths))
B <- getBeta(GRset.funnorm)

data_set_list[["KAT6A"]] <- list( data_set_name = "KAT6A", g1 = g1, g2 = g2, locs = locs, B = B,
                                g1g2_labels = data.frame(g1="Control",g2="KAT6A"))




method_set_list <- list(
  dmrscaler_1 = list(method="dmrscaler",
                     function_call=" DMRscaler::dmrscaler(locs = locs, locs_pval_cutoff = 0.01, region_signif_method = \"bon\", region_signif_cutoff = 0.01, window_type = \"k_near\", window_sizes = c(2,4,8,16,32,64), output_type = \"comp\") " ),
  dmrscaler_2 = list(method="dmrscaler",
                     function_call=" DMRscaler::dmrscaler(locs = locs, locs_pval_cutoff = 0.01, region_signif_method = \"ben\", region_signif_cutoff = 0.01, window_type = \"k_near\", window_sizes = c(2,4,8,16,32,64), output_type = \"comp\") " ),
  bumphunter_1 = list(method="bumphunter",
                      function_call="bumphunter(B,as.matrix(design),chr = locs$chr, pos=locs$pos, cutoff=0.1, maxGap=1e3, B=250, smoothFunction=loessByCluster )"),
  bumphunter_2 = list(method="bumphunter",
                      function_call="bumphunter(B,as.matrix(design),chr = locs$chr, pos=locs$pos, cutoff=0.1, maxGap=1e6, B=250, smoothFunction=loessByCluster )"),
  dmrcate_1 = list(method="dmrcate",
                   function_call="dmrcate(myannotation, lambda=1e3, C=2)"),
  dmrcate_2 = list(method="dmrcate",
                   function_call="dmrcate(myannotation, lambda=1e6, C=2000)"),
  combp_1 = list(method="combp",
                 function_call="system(\"comb-p pipeline -c 4 --dist 1000 --step 100 --seed 1e-3 --region-filter-p 0.1 -p " ),
  combp_2 = list(method="combp",
                 function_call="system(\"comb-p pipeline -c 4 --dist 1000000 --step 5000 --seed 1e-3 --region-filter-p 0.1 -p " )
  combp_3 = list(method="combp",
                 function_call="system(\"comb-p pipeline -c 4 --dist 1000000 --step 100000 --seed 1e-3 --region-filter-p 0.1 -p " )
)


## write Rdata objects to load into each real_data_individual_run.R script
filename <- paste("real_data_setup.Rdata")
save.image(file=filename)
## write data_set_table
data_set_table <- names(data_set_list)
filename <- paste("data_set_table.csv")
write.table(data_set_table, filename, row.names = F, col.names=F, sep=",")

## write method_table
method_table <- names(method_set_list)
filename <- paste("method_table.csv")
write.table(method_table, filename, row.names = F, col.names=F, sep=",")
