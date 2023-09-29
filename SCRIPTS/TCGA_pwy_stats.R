rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel" # pathways preselected as positive controls

type = "gs"
#type = "ph"

pval_thr = c(0.05,0.01,0.001)
n_thr = length(pval_thr)

weight_label = c("cl",paste0("p",c("1","1.5","2")))
n_weight = length(weight_label)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

pval_mat = vector("list",n_proj)
for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_",collection,"_",type,"_pval.txt"))
    data = as.matrix(read.table(infile,header=T,sep="\t"))
    data = data[,-1]
    pval_mat[[i_proj]] = matrix(as.numeric(data),ncol=ncol(data))
}

for (i_thr in 1:n_thr) {
    proj_mat = matrix(rep(NA,n_proj*n_weight),ncol=n_weight)
    for (i_proj in 1:n_proj) {
        proj_mat[i_proj,] = apply(pval_mat[[i_proj]]<pval_thr[i_thr],2,sum)
    }
    outfile = file.path(PROJECT_DIR,"TCGA",paste0("TCGA_",collection,"_",type,"_thr",pval_thr[i_thr],"_pwy_stats.txt"))
    output = rbind(c("Project",paste(type,weight_label,sep="_")),cbind(proj,proj_mat))
    write(t(output),ncol=ncol(output),sep="\t",file=outfile)
}
