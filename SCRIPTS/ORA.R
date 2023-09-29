library(fgsea)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

#collection = "pwy.presel" # pathways preselected as positive controls
collection = "pos.ctrl"
#collection = "neg.ctrl.unif" # randomized negative control pathways with uniform distribution of pathway sizes between size_min=15 and size_max=500

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

method = c("bonferroni","BH")
method_label = c("Bonf","BH")
n_target = length(method)

dir = c("TOP","BOT")
n_dir = length(dir)

size_max = 500
size_min = 15

for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_ranks.rnk"))
    ranked_genes = as.matrix(read.table(infile,sep="\t",header=F))
    remove = duplicated(ranked_genes)
    ranked_genes = ranked_genes[!remove,]
    sign = sign(as.numeric(ranked_genes[,2]))
    pv = 10**-(abs(as.numeric(ranked_genes[,2])))
    BKG = ranked_genes[,1]
    BKG_size = length(BKG)
    
    if (collection=="pwy.presel") {
        gmt_file = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_presel.gmt"))
    } else if (collection=="pos.ctrl") {
        gmt_file = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_pos.ctrl.gmt"))
    } else if (collection=="neg.ctrl.unif") {
        gmt_file = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_pwy_neg.ctrl.unif.gmt"))
    }
    
    ref = fgsea::gmtPathways(gmt_file)
    pwy = names(ref)
    n_pwy = length(pwy)
    pwy_size = rep(NA,n_pwy)
    for (i_pwy in 1:n_pwy) {
        sel = ref[[i_pwy]]%in%BKG
        ref[[i_pwy]] = ref[[i_pwy]][sel]
        pwy_size[i_pwy] = length(ref[[i_pwy]])
    }
    remove = (pwy_size<size_min)|(pwy_size>size_max)
    if (sum(remove)>0) {
        ref = ref[-which(remove)]
        pwy_size = pwy_size[-which(remove)]
        pwy = names(ref)
        n_pwy = length(pwy)
    }
    
    for (i_target in 1:n_target) {
        pval = matrix(rep(NA,n_pwy*n_dir),ncol=n_dir)
        colnames(pval) = dir
        for (i_dir in 1:n_dir) {
            if (dir[i_dir]=="TOP") {
                sel = (sign>0)&(p.adjust(pv,method=method[i_target])<0.05)
                target = BKG[sel]
            } else if (dir[i_dir]=="BOT") {
                sel = (sign<0)&(p.adjust(pv,method=method[i_target])<0.05)
                target = BKG[sel]
            }
            target_size = length(target)
            for (i_pwy in 1:n_pwy) {
                overlap = sum(target%in%ref[[i_pwy]])
                pval[i_pwy,i_dir] = phyper(overlap-1,pwy_size[i_pwy],BKG_size-pwy_size[i_pwy],target_size,lower.tail=F)
            }
        }
        pval_merged = apply(pval,1,min)
        signif_merged = rep(0,n_pwy)
        signif_merged[pval[,"TOP"]<0.05] = 1
        signif_merged[pval[,"BOT"]<0.05] = -1
        signif_merged[(pval[,"TOP"]<0.05)&(pval[,"BOT"]<0.05)] = 2
        
        sel = p.adjust(pv,method=method[i_target])<0.05
        target = BKG[sel]
        target_size = length(target)
        pval_abs = rep(0,n_pwy)
        for (i_pwy in 1:n_pwy) {
            overlap = sum(target%in%ref[[i_pwy]])
            pval_abs[i_pwy] = phyper(overlap-1,pwy_size[i_pwy],BKG_size-pwy_size[i_pwy],target_size,lower.tail=F)
        }
        
        res = cbind(pval,pval_merged,signif_merged,pval_abs)
        output = rbind(c("Pathway",paste0("pval_TOP_",method_label[i_target]),paste0("pval_BOT_",method_label[i_target]),"pval_merged","dir_merged",paste0("pval_abs_",method_label[i_target])),cbind(pwy,res))
        outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_",collection,"_ORA_",method_label[i_target],".txt"))
        write(t(output),ncol=ncol(output),sep="\t",file=outfile)
    }
}
