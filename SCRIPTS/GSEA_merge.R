rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel" # pathways preselected as positive controls
# collection = "neg.ctrl.unif" # randomized negative control pathways with uniform distribution of pathway sizes between size_min=15 and size_max=500

type = "gs"

weight = c("cl",paste0("p",c("1","1.5","2")))
n_weight = length(weight)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
#infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

for (i_proj in 1:n_proj) {
    output_dir = file.path(PROJECT_DIR,proj[i_proj])
    for (i_weight in 1:n_weight) {
        run_label = paste0(proj[i_proj],"_",collection,"_",type,"_",weight[i_weight])
        gsea_dir = list.files(path=output_dir,pattern=paste0(run_label,".GseaPreranked"))
        # I'm assuming only one directory (i.e. GSEA run) matches the pattern
        gsea_files = list.files(path=file.path(output_dir,gsea_dir),pattern="gsea_report_for_")
        gsea_files = gsea_files[grep(gsea_files,pattern="\\.tsv")]
        gsea_res = NULL
        for (i in 1:length(gsea_files)){
            gsea_res_add = as.matrix(read.table(file.path(output_dir,gsea_dir,gsea_files[i]),header=T,sep="\t",quote="\"",stringsAsFactors=F))
            gsea_res = rbind(gsea_res,gsea_res_add)
        }
        if (i_weight==1) {
            pwy = gsea_res[,"NAME"]
            n_pwy = length(pwy)
            dir_mat = matrix(rep(0,n_pwy*n_weight),ncol=n_weight)
            rownames(dir_mat) = pwy
            colnames(dir_mat) = paste(type,weight,sep="_")
            pval_mat = matrix(rep(0,n_pwy*n_weight),ncol=n_weight)
            rownames(pval_mat) = pwy
            colnames(pval_mat) = paste(type,weight,sep="_")
        }
        dir = sign(as.numeric(gsea_res[,"ES"]))
        pval = as.numeric(gsea_res[,"NOM.p.val"])
        o = match(pwy,gsea_res[,"NAME"])
        dir_mat[,i_weight] = dir[o]
        pval_mat[,i_weight] = pval[o]
    }
    # we remove NA's in p-values (rare occurrence).
    index = which(is.na(pval_mat))
    if (length(index)>0) {
        pval_mat[index] = 1
    }
    
    outfile = file.path(output_dir,paste0(proj[i_proj],"_",collection,"_",type,"_dir.txt"))
    output = rbind(c("Pathway",colnames(dir_mat)),cbind(pwy,dir_mat))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
    
    outfile = file.path(output_dir,paste0(proj[i_proj],"_",collection,"_",type,"_pval.txt"))
    output = rbind(c("Pathway",colnames(pval_mat)),cbind(pwy,pval_mat))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
