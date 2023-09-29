library(fgsea)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub" # replace with your local path

infile = file.path(PROJECT_DIR,"RESULTS","TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)
search_terms = vector("list",n_proj)

# manually curated search terms for each cancer type (in the same order as they appear in this file - i.e. alphabetically ordered)
# WARNING: this is hard-coded!
search_terms[[1]] = c("bladder_cancer","bladder_neoplasm")
search_terms[[2]] = c("breast_cancer","breast_carcinoma","neoplasm_of_the_breast")
search_terms[[3]] = c("NEOPLASM_OF_THE_COLON","COLON_CANCER",
"ADENOCARCINOMA_OF_THE_COLON","COLON_AND_RECTAL_CANCER","COLORECTAL_CANCER")
search_terms[[4]] = c("ESOPHAGEAL_CARCINOMA","ESOPHAGEAL_CARCINOGENESIS","ESOPHAGEAL_CARCINOMA",
"ESOPHAGEAL_NEOPLASM","ESOPHAGUS_CANCER")
search_terms[[5]] = c("HEAD_AND_NECK_CANCER","NEOPLASM_OF_HEAD_AND_NECK","HEAD_AND_NECK_SQUAMOUS_CELL_CARCINOMA")
search_terms[[6]] = c("KIDNEY_CANCER","RENAL_CELL_CARCINOMA")
search_terms[[7]] = search_terms[[6]]
search_terms[[8]] = search_terms[[6]]
search_terms[[9]] = c("LIVER_CANCER","HEPATOCELLULAR_CARCINOMA","NEOPLASM_OF_THE_LIVER")
search_terms[[10]] = c("LUNG_CANCER","LUNG_CARCINOMA")
search_terms[[11]] = search_terms[[10]]
search_terms[[12]] = c("PROSTATE_CANCER","PROSTATE_CARCINOGENESIS","PROSTATE_NEOPLASM")
search_terms[[13]] = c("STOMACH_CANCER","GASTRIC_CANCER","NEOPLASM_OF_THE_STOMACH")
search_terms[[14]] = c("_THYROID_CARCINOMA","_THYROID_CANCER","NEOPLASM_OF_THE_THYROID")
search_terms[[15]] = c("ENDOMETRIAL_CANCER")

infile = file.path(PROJECT_DIR,"DATA","msigdb.v2023.1.Hs.symbols.gmt")
ref_pathways = fgsea::gmtPathways(infile)
pwy = names(ref_pathways)

for (i_proj in 1:n_proj) {
    pwy_match = grep(paste0(search_terms[[i_proj]],collapse="|"),pwy,ignore.case=T)
    outfile = file.path(PROJECT_DIR,"RESULTS",proj[i_proj],paste0(proj[i_proj],"_pwy_presel.gmt"))
    writeGmtPathways(ref_pathways[pwy_match],outfile)
}
