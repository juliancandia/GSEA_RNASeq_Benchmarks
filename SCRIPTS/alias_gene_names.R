library(fgsea)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub" # replace with your local path

path_data = "~/DATA/GDC_GEQ_230719" # replace with your local path (location of TCGA files downloaded from the NCI GDC portal - see README file for further details)

infile = file.path(PROJECT_DIR,"RESULTS","TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]

# we open one sample (any sample is OK) - we just want the list of genes
# Note that expression matrices for all samples contain the same set of genes (harmonized pipeline from NCI GDC).
i_proj=1
infile = file.path(PROJECT_DIR,"RESULTS",proj[i_proj],paste0(proj[i_proj],"_sample.txt"))
sample = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
i_sample=1
infile = file.path(path_data,sample[i_sample,"File.ID"],sample[i_sample,"File.Name"])
data = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
gene_metadata = data[-(1:4),1:3]
gene_sel = gene_metadata[,"gene_type"]=="protein_coding"
gene_metadata = gene_metadata[gene_sel,1:2]
colnames(gene_metadata) = c("ENSG.ver","name")
ENSG = matrix(unlist(strsplit(gene_metadata[,"ENSG.ver"],"\\.")),byrow=T,ncol=2)[,1]
gene_metadata = cbind(ENSG,gene_metadata)

# aliases to match pathway names
infile = file.path(PROJECT_DIR,"DATA","genenames_map.txt")
map = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t",quote="",comment.char=""))
infile = file.path(PROJECT_DIR,"DATA","msigdb.v2023.1.Hs.symbols.gmt")
ref_pathways = fgsea::gmtPathways(infile) # 33,591 pathways
ref_genes = unique(unlist(ref_pathways)) # 42228 genes

alias = gene_metadata[,"name"]
sel_missing = !alias%in%ref_genes
gene_missing = gene_metadata[sel_missing,c("ENSG","name")]
n_missing_pre = nrow(gene_missing)
for (i in 1:nrow(gene_missing)) {
    symbol = map[map[,"ENSG"]==gene_missing[i,"ENSG"],"Symbol"]
    sel = symbol%in%ref_genes
    if (sum(sel)>0) {
        gene_missing[i,"name"] = symbol[sel][1]
    }
}
alias[sel_missing] = gene_missing[,"name"]
n_missing_post = sum(!gene_missing[,"name"]%in%ref_genes)
cat(paste0("Number of gene aliases found:"),n_missing_pre-n_missing_post,"\n")
in_msigdb = rep("Y",nrow(gene_metadata))
in_msigdb[!alias%in%ref_genes] = "N"
gene_metadata = cbind(gene_metadata,alias,in_msigdb)

# we save the annotations
outfile = file.path(PROJECT_DIR,"RESULTS","TCGA","gene_metadata.txt")
output = rbind(colnames(gene_metadata),gene_metadata)
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
