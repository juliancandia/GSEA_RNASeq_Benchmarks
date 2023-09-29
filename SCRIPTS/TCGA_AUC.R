rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

#collection = c("pwy.presel","neg.ctrl.unif") # (TRUE,NEG) collections to compare
collection = c("pos.ctrl","neg.ctrl.unif") # (TRUE,NEG) collections to compare

weight = c("cl",paste0("p",c("1","1.5","2")))
n_weight = length(weight)
type = c("gs","ph")
type_label = c("gene-set permutations","phenotype permutations")
n_type = length(type)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)
proj_label = matrix(unlist(strsplit(proj,"_")),byrow=T,ncol=2)[,2] # short version

for (i_type in 1:n_type) {
    auc_low = matrix(rep(NA,n_proj*n_weight),ncol=n_weight)
    auc = matrix(rep(NA,n_proj*n_weight),ncol=n_weight)
    auc_high = matrix(rep(NA,n_proj*n_weight),ncol=n_weight)
    for (i_proj in 1:n_proj) {
        proj_dir = file.path(PROJECT_DIR,proj[i_proj])
        infile =
        file.path(proj_dir,paste0(proj[i_proj],"_",type[i_type],"_AUC_",paste0(collection,collapse="_"),".txt"))
        data = as.matrix(read.table(infile,sep="\t",header=T,stringsAsFactors=F))
        auc_low[i_proj,] = as.numeric(data[,"auc_low"])
        auc[i_proj,] = as.numeric(data[,"auc"])
        auc_high[i_proj,] = as.numeric(data[,"auc_high"])
    }
    outfile = file.path(PROJECT_DIR,"TCGA",paste0("TCGA_",type[i_type],"_AUC_",paste0(collection,collapse="_"),".pdf"))
    pdf(outfile,width=7,height=5)
    x = 1:n_proj
    range_x = range(x)
    range_y = c(0,1)
    title = paste0("GSEA (",type_label[i_type],"): ",collection[1]," vs ",collection[2])
    plot(range_x,range_y,xaxt="n",type="n",xlab="",ylab="ROC Area Under Curve (with 95% CI)",cex.lab=1,main=title)
    col = c(1:n_weight)
    pch = c(15,16,17,18)
    offset = 0.15
    offset_seq = seq(-offset,offset,length.out=n_weight)
    for (i_weight in 1:n_weight) {
        points(x+offset_seq[i_weight],auc[,i_weight],cex=1,col=col[i_weight],pch=pch[i_weight],lty=1,type="b")
        arrows(x+offset_seq[i_weight],auc_low[,i_weight],x+offset_seq[i_weight],
        auc_high[,i_weight],col=col[i_weight],length=0.025,angle=90,code=3)
    }
    axis(side=1,at=x,labels=proj_label,las=2,cex.axis=1)
    abline(h=0.5,col="black",lty=3)
    legend("bottomleft",weight,pch=pch,lwd=1,col=col,lty=1,cex=1,pt.cex=1)
    dev.off()
}
