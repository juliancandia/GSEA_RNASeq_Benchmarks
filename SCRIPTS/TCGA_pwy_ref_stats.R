library(fgsea)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

pwy_stats = matrix(rep(NA,n_proj*3),ncol=3)
for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_presel.gmt"))
    ref = fgsea::gmtPathways(infile)
    pwy_stats[i_proj,1] = length(ref)
    
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy.presel_gs_pval.txt"))
    data = as.matrix(read.table(infile,header=T,sep="\t"))
    pwy_stats[i_proj,2] = nrow(data)
    
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_pos.ctrl.gmt"))
    ref = fgsea::gmtPathways(infile)
    pwy_stats[i_proj,3] = length(ref)
}

outfile = file.path(PROJECT_DIR,"TCGA","TCGA_pwy_ref_stats.txt")
output = rbind(c("Project","Preselected","Target","Positive-Control"),cbind(proj,pwy_stats))
write(t(output),ncol=ncol(output),sep="\t",file=outfile)
