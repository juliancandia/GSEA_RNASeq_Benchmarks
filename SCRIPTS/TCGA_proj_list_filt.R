rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel" # pathways preselected as positive controls

type = "gs" # most sensitive permutation method

pval_thr = 0.05 # most lenient

infile = file.path(PROJECT_DIR,"TCGA",paste0("TCGA_",collection,"_",type,"_thr",pval_thr,"_pwy_stats.txt"))
data = as.matrix(read.table(infile,header=T,sep="\t"))
proj = data[,1]
data = data[,-1]
signif_pwys = matrix(as.numeric(data),ncol=ncol(data))
remove = apply(signif_pwys,1,mean)<1
proj = proj[!remove]

outfile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
write(proj,ncol=1,sep="\t",file=outfile)
