rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","sample_metadata.txt")
sample = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
sel = grep("TCGA",sample[,"Project.ID"])
sample = sample[sel,] # 11,274 TCGA samples (downloaded from NCI's GDC).

proj = unique(sample[,"Project.ID"])
proj = proj[order(proj)]
proj_label = paste0("TCGA_",matrix(unlist(strsplit(proj,"-")),ncol=2,byrow=T)[,2])
# NOTE: we replace hyphens by underscores because GSEA can't handle hyphens
n_proj = length(proj)
proj_sel = rep(F,n_proj)
for (i_proj in 1:n_proj) {
    sample_this_proj = sample[sample[,"Project.ID"]==proj[i_proj],]
    sel = sample_this_proj[,"Sample.Type"]%in%c("Solid Tissue Normal","Primary Tumor")
    sample_this_proj = sample_this_proj[sel,]
    sel_N = sample_this_proj[,"Sample.Type"]=="Solid Tissue Normal"
    ID_N = unique(sample_this_proj[sel_N,"Case.ID"])
    sel_T = sample_this_proj[,"Sample.Type"]=="Primary Tumor"
    ID_T = unique(sample_this_proj[sel_T,"Case.ID"])
    # paired analysis
    ID_NT = ID_N[ID_N%in%ID_T]
    if (length(ID_NT)>=10) { # we require at least 10 subjects with paired tumor-normal samples
        sample_this_proj = sample_this_proj[sample_this_proj[,"Case.ID"]%in%ID_NT,]
        output_dir = file.path(PROJECT_DIR,"RESULTS",proj_label[i_proj])
        dir.create(output_dir,showWarnings=F,recursive=T)
        output = rbind(colnames(sample_this_proj),sample_this_proj)
        outfile = file.path(output_dir,paste0(proj_label[i_proj],"_sample.txt"))
        write(t(output),ncol=ncol(output),file=outfile,sep="\t")
        proj_sel[i_proj] = T
    }
}

output_dir = file.path(PROJECT_DIR,"RESULTS","TCGA")
dir.create(output_dir,showWarnings=F)
outfile = file.path(output_dir,"TCGA_proj_list.txt")
proj = proj_label[proj_sel]
write(proj,ncol=1,file=outfile,sep="\t")
