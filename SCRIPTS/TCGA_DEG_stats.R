rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

n_sample_analyzed = rep(NA,n_proj)
n_subj_analyzed = rep(NA,n_proj)
n_gene = rep(NA,n_proj)
n_pwy = rep(NA,n_proj)
res = NULL
for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_ranks.rnk"))
    rank = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))
    dir = sign(as.numeric(rank[,2]))
    pval = 10**-(abs(as.numeric(rank[,2])))
    bh = p.adjust(pval,method="BH")
    bonf = p.adjust(pval,method="bonferroni")
    res = rbind(res,c(length(dir),sum(bh<0.05),sum(bonf<0.05)))
}
output = rbind(c("Project","Genes","BH<0.05","Bonferroni<0.05"),cbind(proj,res))
outfile = file.path(PROJECT_DIR,"TCGA","TCGA_DEG_stats.txt")
write(t(output),ncol=ncol(output),sep="\t",file=outfile)
