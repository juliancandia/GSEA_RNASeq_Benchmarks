library(fgsea)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection_input = "pwy.presel"
collection_output = "pos.ctrl"

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

pval_thr = 0.05

type = c("gs","ph")
n_type = length(type)

for (i_proj in 1:n_proj) {
    pwy_mat = vector("list",n_type)
    for (i_type in 1:n_type) {
        infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_",collection_input,"_",type[i_type],"_pval.txt"))
        pwy_mat[[i_type]] = as.matrix(read.table(infile,header=T,sep="\t"))
        if (i_type==1) {
            pwy = pwy_mat[[i_type]][,"Pathway"]
            pwy_pval = pwy_mat[[i_type]][,-1]
            pwy_pval_merged = matrix(as.numeric(pwy_pval),ncol=ncol(pwy_pval))
        } else {
            o = match(pwy,pwy_mat[[i_type]][,"Pathway"])
            pwy_mat[[i_type]] = pwy_mat[[i_type]][o,]
            cat("Checkpoint passed:",identical(pwy,pwy_mat[[i_type]][,"Pathway"]),"\n")
            pwy_pval = pwy_mat[[i_type]][,-1]
            pwy_pval_merged = cbind(pwy_pval_merged,matrix(as.numeric(pwy_pval),ncol=ncol(pwy_pval)))
        }
    }
    # we filter preselected pathways by positive-control criteria
    pwy_sel = apply(pwy_pval_merged<pval_thr,1,sum)>0
    for (i_type in 1:n_type) {
        outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_",collection_output,"_",type[i_type],"_pval.txt"))
        output = rbind(colnames(pwy_mat[[i_type]]),pwy_mat[[i_type]][pwy_sel,])
        write(t(output),ncol=ncol(output),file=outfile,sep="\t")
    }
    
    # we generate the filtered gmt file
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_presel.gmt"))
    ref = fgsea::gmtPathways(infile)
    ref = ref[names(ref)%in%pwy[pwy_sel]]
    outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_",collection_output,".gmt"))
    writeGmtPathways(ref,outfile)
}
