library(fgsea)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

size_min = 15
size_max = 500
n_rdm = 1000

for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_ranks.rnk"))
    ranked_genes = as.matrix(read.table(infile,sep="\t",header=F))
    BKG = unique(ranked_genes[,1])
    size = round(runif(n=n_rdm,min=size_min,max=size_max))
    pwy_rdm = vector("list",n_rdm)
    for (i_rdm in 1:n_rdm) {
        pwy_rdm[[i_rdm]] = sample(BKG,size[i_rdm])
    }
    names(pwy_rdm) = paste0(proj[i_proj],"_neg.ctrl.unif.",1:n_rdm)
    outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_neg.ctrl.unif.gmt"))
    writeGmtPathways(pwy_rdm,outfile)
}
