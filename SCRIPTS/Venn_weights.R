library(VennDiagram)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel" # pathways preselected as positive controls
#collection = "pos.ctrl"

#pval_thr = c(0.05,0.01,0.001)
pval_thr = c(0.05)
n_thr = length(pval_thr)

weight_label_full = c("classic",paste0("weighted_",paste0("p",c("1","1.5","2"))))
weight_label = c("cl",paste0("p",c("1","1.5","2")))
n_weight = length(weight_label)

type = c("gs","ph")
n_type = length(type)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)
for (i_proj in 1:n_proj) {
    output_dir = file.path(PROJECT_DIR,proj[i_proj])
    pwy_mat = vector("list",n_type)
    for (i_type in 1:n_type) {
        infile = file.path(output_dir,paste0(proj[i_proj],"_",collection,"_",type[i_type],"_pval.txt"))
        pwy_mat[[i_type]] = as.matrix(read.table(infile,header=T,sep="\t"))
    }
    for (i_weight in 1:n_weight) {
        for (i_thr in 1:n_thr) {
            outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_",collection,"_",weight_label[i_weight],"_Venn_thr",pval_thr[i_thr],".tiff"))
            pwy_type = vector("list",n_type)
            for (i_type in 1:n_type) {
                pval = as.numeric(pwy_mat[[i_type]][,paste0(type[i_type],"_",weight_label[i_weight])])
                sel = pval<pval_thr[i_thr]
                pwy_type[[i_type]] = pwy_mat[[i_type]][sel,"Pathway"]
            }
            type_color = rainbow(n_type)
            if (length(unlist(pwy_type))>0) {
                title = paste0(weight_label_full[i_weight]," (pval<",pval_thr[i_thr],")")
                venn.diagram(x=pwy_type,filename=outfile,main=title,category.names=type,fill=type_color,cat.cex=2,cat.fontface="bold",main.cex=2,cex=2,disable.logging=T)
            }
        }
    }
}
