
rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel" # pathways preselected as positive controls
# collection = "neg.ctrl.unif" # randomized negative control pathways with uniform distribution of pathway sizes between size_min=15 and size_max=500

weight_label = c("cl",paste0("p",c("1","1.5","2")))
n_weight = length(weight_label)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
#infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

for (i_proj in 1:n_proj) {
    output_dir = file.path(PROJECT_DIR,proj[i_proj])
    infile = file.path(output_dir,paste0(proj[i_proj],"_",collection,"_gs_pval.txt"))
    pwy_mat_gs = as.matrix(read.table(infile,header=T,sep="\t"))
    pwy = pwy_mat_gs[,"Pathway"]
    n_pwy = length(pwy)
    pwy_mat_ph = matrix(rep(0,n_pwy*n_weight),ncol=n_weight)
    for (i_weight in 1:n_weight) {
        # We extract ES from the GSEA true (non-shuffled) data
        run_label = paste0(proj[i_proj],"_",collection,"_gs_",weight_label[i_weight])
        gsea_dir = list.files(path=output_dir,pattern=paste0("^",run_label,".GseaPreranked"))
        gsea_files = list.files(path=file.path(output_dir,gsea_dir),pattern="gsea_report_for_na")
        gsea_files = gsea_files[grep(gsea_files,pattern="\\.tsv")]
        gsea_res = NULL
        for (i in 1:length(gsea_files)){
            gsea_res_add = as.matrix(read.table(file.path(output_dir,gsea_dir,gsea_files[i]),header=T,sep="\t",quote="\"",stringsAsFactors=F))
            gsea_res = rbind(gsea_res,gsea_res_add)
        }
        gsea_res = gsea_res[match(pwy,gsea_res[,"NAME"]),]
        # We extract ES from the shuffled data
        output_dir_rdm = file.path(output_dir,paste0(proj[i_proj],"_rdm_",collection))
        run_label_rdm = paste0(proj[i_proj],"_",collection,"_rdm_",weight_label[i_weight])
        rdm_dir = list.files(path=output_dir_rdm,pattern=paste0("^",run_label_rdm,".GseaPreranked"))
        n_rdm = length(rdm_dir)
        rdm_ES = matrix(rep(NA,n_pwy*n_rdm),ncol=n_rdm)
        for (j in 1:n_rdm) {
            rdm_files = list.files(path=file.path(output_dir_rdm,rdm_dir[j]),pattern="gsea_report_for_na")
            rdm_files = rdm_files[grep(rdm_files,pattern="\\.tsv")]
            rdm_res = NULL
            for (i in 1:length(rdm_files)){
                rdm_res_add = as.matrix(read.table(file.path(output_dir_rdm,rdm_dir[j],rdm_files[i]),header=T,sep="\t",quote="\"",stringsAsFactors=F))
                rdm_res = rbind(rdm_res,rdm_res_add)
            }
            o = match(pwy,rdm_res[,"NAME"])
            rdm_ES[,j] = as.numeric(rdm_res[o,"ES"])
        }
        pwy_pval = rep(NA,n_pwy)
        for (i_pwy in 1:n_pwy) {
            pwy_pval[i_pwy] = (sum(abs(rdm_ES[i_pwy,])>abs(as.numeric(gsea_res[i_pwy,"ES"])))+1)/(n_rdm+1)
        }
        pwy_mat_ph[,i_weight] = pwy_pval
    }
    outfile = file.path(output_dir,paste0(proj[i_proj],"_",collection,"_ph_pval.txt"))
    output = rbind(c("Pathway",paste("ph",weight_label,sep="_")),cbind(pwy,pwy_mat_ph))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
