library(gplots)
library(RColorBrewer)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel"

pval_thr = 0.05
type = "gs"
output_prefix = "REBC_THCA"
proj = c("TCGA_THCA","REBC_THYR")
n_proj = length(proj)

EES = vector("list",n_proj) # Enrichment Evidence Score
for (i_proj in 1:n_proj) {
    proj_dir = file.path(PROJECT_DIR,proj[i_proj])
    
    infile = file.path(proj_dir,paste0(proj[i_proj],"_",collection,"_",type,"_pval.txt"))
    pwy_pval = as.matrix(read.table(infile,header=T,sep="\t",stringsAsFactors=F))
    pwy = pwy_pval[,"Pathway"]
    n_pwy = length(pwy)
    col_index = grep(type,colnames(pwy_pval))
    pwy_pval_num = pwy_pval[,col_index]
    pwy_pval_num = matrix(as.numeric(pwy_pval_num),ncol=ncol(pwy_pval_num))
    pwy_pval_bin = (pwy_pval_num<pval_thr)*1
    
    infile = file.path(proj_dir,paste0(proj[i_proj],"_",collection,"_",type,"_dir.txt"))
    pwy_dir = as.matrix(read.table(infile,header=T,sep="\t",stringsAsFactors=F))
    pwy_dir_num = pwy_dir[,col_index]
    pwy_dir_num = matrix(as.numeric(pwy_dir_num),ncol=ncol(pwy_dir_num))
    
    EES[[i_proj]] = apply(pwy_pval_bin*pwy_dir_num,1,sum)
    names(EES[[i_proj]]) = pwy
}

# NOTE: we can't assume the pathways are exactly the same, because they are filtered by size, and that is affected by the set of measured genes.
pwy = unique(c(names(EES[[1]]),names(EES[[2]])))
EES_merged = matrix(rep(NA,2*length(pwy)),ncol=2)
for (i_proj in 1:n_proj) {
    o = match(pwy,names(EES[[i_proj]]))
    EES_merged[,i_proj] = EES[[i_proj]][o]
}
tmp = EES_merged
index = which(is.na(tmp),arr.ind=T)
if (length(index)>0) {
    tmp[index] = 0
}
o = order(-apply(abs(tmp),1,sum))
o2 = order(apply(is.na(EES_merged[o,]),1,sum))
o = o[o2]
output = rbind(c("Pathway",proj),cbind(pwy,EES_merged)[o,])
outfile = file.path(PROJECT_DIR,"TCGA",paste0(output_prefix,"_EES.txt"))
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

sel = names(EES[[1]])%in%names(EES[[2]])
EES[[1]] = EES[[1]][sel]
sel = names(EES[[2]])%in%names(EES[[1]])
EES[[2]] = EES[[2]][sel]
pwy = names(EES[[1]])
for (i_proj in 2:n_proj) {
    o = match(pwy,names(EES[[i_proj]]))
    cat("Checkpoint passed:",identical(pwy,names(EES[[i_proj]])[o]),"\n")
    EES[[i_proj]] = EES[[i_proj]][o]
}

score = seq(-4,4,1)
n_score = length(score)
dt = matrix(rep(0,n_score**2),ncol=n_score)
for (i_score in 1:n_score) {
    for (j_score in 1:n_score) {
        dt[i_score,j_score] = sum((EES[[1]]==score[i_score])&(EES[[2]]==score[j_score]))
    }
}
rownames(dt) = score
colnames(dt) = score

# Balloonplot (full EES matrix)
outfile = file.path(PROJECT_DIR,"TCGA",paste0(output_prefix,"_comp_full.pdf"))
pdf(outfile,width=5,height=5)
dt = as.table(dt)
n = nrow(dt)
dist_mat = matrix(rep(0,n**2),nrow=n)
for (i in 1:(n-1)) {
    for (j in (i+1):n) {
        dist_mat[i,j] = j-i
        dist_mat[j,i] = j-i
    }
}
col = brewer.pal(7,"Blues")[-1]
col_mat = dist_mat
for (i in 1:n) {
    for (j in 1:n) {
        col_mat[i,j] = col[dist_mat[i,j]+1]
    }
}
balloonplot(t(dt),main="",xlab="",ylab="",label=T,show.margins=F,dotcolor=as.character(col_mat),dotsize=5)
dev.off()

res = fisher.test(dt,simulate.p.value=TRUE)
cat(paste0("Full matrix: Fisher's Exact Test p-value (simulated)=",signif(res$p.value,digits=2)),"\n")

score_label = c("Non-tumor","Undet","Tumor")
score_min = c(-4,-2,3)
score_max = c(-3,2,4)
n_score = length(score_label)
dt = matrix(rep(0,n_score**2),ncol=n_score)
for (i_score in 1:n_score) {
    for (j_score in 1:n_score) {
        dt[i_score,j_score] = sum((EES[[1]]>=score_min[i_score])&(EES[[1]]<=score_max[i_score])&(EES[[2]]>=score_min[j_score])&(EES[[2]]<=score_max[j_score]))
    }
}
rownames(dt) = score_label
colnames(dt) = score_label

# Balloonplot (grouped EES matrix)
outfile = file.path(PROJECT_DIR,"TCGA",paste0(output_prefix,"_comp_grouped.pdf"))
pdf(outfile,width=5,height=5)
dt = as.table(dt)
n = nrow(dt)
dist_mat = matrix(rep(0,n**2),nrow=n)
for (i in 1:(n-1)) {
    for (j in (i+1):n) {
        dist_mat[i,j] = j-i
        dist_mat[j,i] = j-i
    }
}
col = brewer.pal(4,"Blues")[-1]
col_mat = dist_mat
for (i in 1:n) {
    for (j in 1:n) {
        col_mat[i,j] = col[dist_mat[i,j]+1]
    }
}
balloonplot(t(dt),main="",xlab="",ylab="",label=T,show.margins=F,dotcolor=as.character(col_mat),dotsize=12)
dev.off()

res = fisher.test(dt)
cat(paste0("Grouped-EES matrix: Fisher's Exact Test p-value=",signif(res$p.value,digits=2)),"\n")
