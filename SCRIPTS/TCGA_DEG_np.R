library(edgeR)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

path_data = "~/DATA/GDC_GEQ_230719" # replace with your local path (location of TCGA files downloaded from the NCI GDC portal - see README file for further details)

infile = file.path(PROJECT_DIR,"TCGA","gene_metadata.txt")
gene_metadata = as.matrix(read.table(infile,header=T,sep="\t"))

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

# we check to see if we can remove FFPE samples
for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_sample.txt"))
    sample = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
    ffpe_index = which(sample[,"is_ffpe"]=="True")
    notffpe_index = which(sample[,"is_ffpe"]=="False")
    n_index = length(ffpe_index)
    if (n_index>0) {
        for (i_index in 1:n_index) {
            subj = sample[ffpe_index[i_index],"Case.ID"]
            type = sample[ffpe_index[i_index],"Sample.Type"]
            notffpe_exists = sum((sample[notffpe_index,"Case.ID"]==subj)&(sample[notffpe_index,"Sample.Type"]==type))>0
            if (notffpe_exists) {
                cat("Proj=",proj[i_proj],"/ Subj=",subj,"/ OK","\n")
            } else {
                cat("Proj=",proj[i_proj],"/ Subj=",subj,"/ FFPE only","\n")
            }
        }
    }
}
# the results are that FFPE samples are redundant and can be removed from the analysis

for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_sample.txt"))
    sample = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
    sel = !toupper(sample[,"is_ffpe"])=="TRUE"
    sample = sample[sel,] # we remove FFPE samples
    n_sample = nrow(sample)
    for (i_sample in 1:n_sample) {
        infile = file.path(path_data,sample[i_sample,"File.ID"],sample[i_sample,"File.Name"])
        data = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
        if (i_sample==1) {
            gene_annot = data[-(1:4),1:3]
            gene_sel = gene_annot[,"gene_type"]=="protein_coding"
            gene_annot = gene_annot[gene_sel,1:2]
            n_gene = nrow(gene_annot)
            expr = matrix(rep(NA,n_gene*n_sample),ncol=n_sample)
        }
        expr[,i_sample] = as.numeric(data[-(1:4),"tpm_unstranded"])[gene_sel]
    }
    keep = rowSums(expr>1) >= round(0.4*n_sample) # arbitrary cutoffs to remove low expressed genes (this is similar to SupplProt4 in Nat Protocols Reimand 2019). Note that the TPM count matrix is normalized to 1 million per sample (across all genes) and it is given as a continuous range of values >=0.
    gene_annot = gene_annot[keep,]
    gene_annot = gene_metadata[match(gene_annot[,"gene_id"],gene_metadata[,"ENSG.ver"]),]
    expr = log2(expr[keep,]+1)
    subj = sample[,"Case.ID"]
    subj_u = unique(subj)
    n_subj = length(subj_u)
    for (i_subj in 1:n_subj) {
        subj[subj==subj_u[i_subj]] = paste0("subj_",i_subj)
    }
    type = sample[,"Sample.Type"]
    type[type=="Solid Tissue Normal"] = "N"
    type[type=="Primary Tumor"] = "T"
    mm = model.matrix(~ 0 + type) # unpaired analysis
    fit <- lmFit(expr,mm)
    contr <- makeContrasts(typeT - typeN, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp,trend=TRUE)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    rank_metric = sign(top.table[,"logFC"])*(-log10(top.table[,"P.Value"]))
    o = order(-rank_metric)
    output = cbind(gene_annot[as.numeric(rownames(top.table)),"alias"],rank_metric)[o,]
    outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_ranks_np.rnk"))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
