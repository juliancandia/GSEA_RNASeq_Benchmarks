library(edgeR)
library(foreach)
library(doParallel)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

n_rdm = 1000

for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_genes.txt"))
    gene_annot = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_sample_annot.txt"))
    sample_annot = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_expr.txt"))
    expr = as.matrix(read.table(infile,header=F,sep="\t"))
    type = sample_annot[,"type"]
    # parallel randomization
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    output_dir_rdm = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_rdm_ranks_np"))
    dir.create(output_dir_rdm,showWarnings=F)
    data = foreach (i_rdm=1:n_rdm,.packages="edgeR",.combine=rbind) %dopar% {
        try({
            type_rdm = sample(type) # unpaired randomization
            mm = model.matrix(~ 0 + type_rdm) # unpaired analysis
            fit <- lmFit(expr,mm)
            contr <- makeContrasts(type_rdmT - type_rdmN, levels = colnames(coef(fit)))
            tmp <- contrasts.fit(fit, contr)
            tmp <- eBayes(tmp,trend=TRUE)
            top.table <- topTable(tmp, sort.by = "P", n = Inf)
            rank_metric = sign(top.table[,"logFC"])*(-log10(top.table[,"P.Value"]))
            o = order(-rank_metric)
            output = cbind(gene_annot[as.numeric(rownames(top.table)),"alias"],rank_metric)[o,]
            outfile = file.path(output_dir_rdm,paste0("rdm",i_rdm,"_ranks_np.rnk"))
            write(t(output),ncol=ncol(output),file=outfile,sep="\t")
        })
    }
    stopCluster(cl)
}
