library(edgeR)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

make_plot = T

header = c("Measured","Bonf<0.05 (FC<0)","Bonf<0.05 (FC>0)","BH<0.05 (FC<0)","BH<0.05 (FC>0)")
gene_stats = matrix(rep(NA,5*n_proj),ncol=5)
for (i_proj in 1:n_proj) {
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_genes.txt"))
    gene_annot = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_sample_annot.txt"))
    sample_annot = as.matrix(read.table(infile,header=T,sep="\t",quote="\"",stringsAsFactors=F))
    infile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_expr.txt"))
    expr = as.matrix(read.table(infile,header=F,sep="\t"))
    subj = sample_annot[,"subj"]
    type = sample_annot[,"type"]
    mm = model.matrix(~ 0 + type + subj) # paired analysis
    fit <- lmFit(expr,mm)
    contr <- makeContrasts(typeT - typeN, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp,trend=TRUE)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    rownames(top.table) = gene_annot[as.numeric(rownames(top.table)),"alias"]
    
    x0 = top.table[,"logFC"]
    y0 = -log10(top.table[,"P.Value"])
    bonf = p.adjust(top.table[,"P.Value"],method="bonferroni")
    bh = p.adjust(top.table[,"P.Value"],method="BH")
    pval_thres = 0.05
    if (make_plot) {
        outfile = file.path(PROJECT_DIR,proj[i_proj],paste0(proj[i_proj],"_Volcano.pdf"))
        pdf(outfile,width=5,height=5)
        x_range = range(x0)
        y_range = range(y0)
        plot(x_range,y_range,type="n",xlab="logFC (Tumor/Non-Tumor)",ylab="-log10(p-value)",cex.lab=1.,main=paste(proj[i_proj]," (paired DGE analysis)"))
        transparency = 80
        up_col = "firebrick3"
        down_col = "navy"
        ns_col = "lightgrey"
        print_labels = F
        n_labels = 10
        sel = bonf>pval_thres
        col_par = as.numeric(col2rgb(ns_col))
        lines(x0[sel],y0[sel],type="p",pch=16,col=rgb(col_par[1],col_par[2],col_par[3],
        transparency,maxColorValue=255))
        sel = (bonf<=pval_thres)&(x0>0)
        col_par = as.numeric(col2rgb(up_col))
        lines(x0[sel],y0[sel],type="p",pch=16,col=rgb(col_par[1],col_par[2],col_par[3],
        transparency,maxColorValue=255))
        sel = (bonf<=pval_thres)&(x0<0)
        col_par = as.numeric(col2rgb(down_col))
        lines(x0[sel],y0[sel],type="p",pch=16,col=rgb(col_par[1],col_par[2],col_par[3],
        transparency,maxColorValue=255))
        border_top = -log10(max(top.table[bonf<=pval_thres,"P.Value"]))
        border_bot = -log10(min(top.table[bonf>pval_thres,"P.Value"]))
        border = mean(c(border_top,border_bot))
        abline(h=border,lty=2)
        abline(v=0,lty=2)
        if (print_labels) {
            o = order(-y0)
            # top positive
            index = o[x0[o]>0][1:n_labels]
            textxy(x0[index],y0[index],rownames(top.table)[index],cex = 0.5, offset = 0.65)
            # top negative
            index = o[x0[o]<0][1:n_labels]
            textxy(x0[index],y0[index],rownames(top.table)[index],cex = 0.5, offset = 0.65)
        }
        dev.off()
    }
    gene_stats[i_proj,1] = nrow(top.table)
    gene_stats[i_proj,2] = sum((bonf<=pval_thres)&(x0<0))
    gene_stats[i_proj,3] = sum((bonf<=pval_thres)&(x0>0))
    gene_stats[i_proj,4] = sum((bh<=pval_thres)&(x0<0))
    gene_stats[i_proj,5] = sum((bh<=pval_thres)&(x0>0))
}

output = rbind(c("Project",header),cbind(proj,gene_stats))
outfile = file.path(PROJECT_DIR,"TCGA","TCGA_DEG_stats_v2.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
