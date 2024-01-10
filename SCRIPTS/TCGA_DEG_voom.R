library(edgeR)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK" # replace with your local path

infile = file.path(PROJECT_DIR,"RESULTS","TCGA","gene_metadata.txt")
gene_metadata = as.matrix(read.table(infile,header=T,sep="\t"))

proj = "TCGA_BRCA"
infile = file.path(PROJECT_DIR,"RESULTS",proj,paste0(proj,"_sample_annot.txt"))
sample = as.matrix(read.table(infile,header=T,sep="\t"))
subj = sample[,"subj"]
type = sample[,"type"]
n_sample = nrow(sample)
for (i_sample in 1:n_sample) {
    infile = file.path(PROJECT_DIR,"REVISION","DATA",sample[i_sample,"File.Name"])
    data = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
    if (i_sample==1) {
        gene_annot = data[-(1:4),1:3]
        gene_sel = gene_annot[,"gene_type"]=="protein_coding"
        gene_annot = gene_annot[gene_sel,1:2]
        n_gene = nrow(gene_annot)
        expr = matrix(rep(NA,n_gene*n_sample),ncol=n_sample)
    }
    expr[,i_sample] = as.numeric(data[-(1:4),"unstranded"])[gene_sel]
}
counts = data.frame(expr)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
# using function from edgeR
keep = filterByExpr(d0, group=type)
d = d0[keep,]
gene_metadata = gene_metadata[keep,]
mm = model.matrix(~ 0 + type + subj)
y0 <- voom(d, mm, plot=F, normalize.method="quantile")
fit <- lmFit(y0,mm)
contr <- makeContrasts(typeT - typeN, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp,trend=TRUE)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
rank_metric = sign(top.table[,"logFC"])*(-log10(top.table[,"P.Value"]))
o = order(-rank_metric)
output = cbind(gene_annot[as.numeric(rownames(top.table)),"gene_name"],rank_metric)[o,]
outfile = file.path(PROJECT_DIR,"REVISION","RESULTS",paste0(proj,"_ranks.rnk"))
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
