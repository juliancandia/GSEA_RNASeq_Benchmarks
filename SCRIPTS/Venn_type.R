library(VennDiagram)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel" # pathways preselected as positive controls
#collection = "pos.ctrl"

#type = "gs"
#type_label = "gene-set permutation"
type = "ph"
type_label = "phenotype permutation"

pval_thr = c(0.05,0.01,0.001)
n_thr = length(pval_thr)

weight_label = c("cl",paste0("p",c("1","1.5","2")))
n_weight = length(weight_label)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)
for (i_proj in 1:n_proj) {
    output_dir = file.path(PROJECT_DIR,proj[i_proj])
    infile = file.path(output_dir,paste0(proj[i_proj],"_",collection,"_",type,"_pval.txt"))
    pwy_mat = as.matrix(read.table(infile,header=T,sep="\t"))
    for (i_thr in 1:n_thr) {
        outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_",collection,"_",type,"_Venn_thr",pval_thr[i_thr],".tiff"))
        pwy_weight = vector("list",n_weight)
        for (i_weight in 1:n_weight) {
            pval = as.numeric(pwy_mat[,paste0(type,"_",weight_label[i_weight])])
            sel = pval<pval_thr[i_thr]
            pwy_weight[[i_weight]] = pwy_mat[sel,"Pathway"]
        }
        weight_color = rainbow(n_weight)
        if (length(unlist(pwy_weight))>0) {
            title = paste0(type_label," (pval<",pval_thr[i_thr],")")
            venn.diagram(x=pwy_weight,filename=outfile,main=title,category.names=weight_label,fill=weight_color,cat.cex=2,cat.fontface="bold",main.cex=2,cex=1.7,disable.logging=T,resolution = 500)
        }
    }
}
