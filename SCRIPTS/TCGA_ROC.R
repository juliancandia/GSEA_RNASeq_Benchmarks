library(pROC)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

type = "gs"
#type = "ph"

#collection = c("pwy.presel","neg.ctrl.unif") # (TRUE,NEG) collections to compare
collection = c("pos.ctrl","neg.ctrl.unif") # (TRUE,NEG) collections to compare
n_coll = 2

weight = c("cl",paste0("p",c("1","1.5","2")))
n_weight = length(weight)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

for (i_proj in 1:n_proj) {
    output_dir = file.path(PROJECT_DIR,proj[i_proj])
    pval_mat = vector("list",n_coll)
    for (i_coll in 1:n_coll) {
        infile = file.path(output_dir,paste0(proj[i_proj],"_",collection[i_coll],"_",type,"_pval.txt"))
        pval_mat[[i_coll]] = as.matrix(read.table(infile,header=T,stringsAsFactors=F,sep="\t"))
    }
    
    # we generate ROC plots / AUC results
    outfile = file.path(output_dir,paste0(proj[i_proj],"_",type,"_ROC_",paste0(collection,collapse="_"),".pdf"))
    pdf(outfile,width=7,height=5)
    true = c(rep(1,nrow(pval_mat[[1]])),rep(0,nrow(pval_mat[[2]])))
    auc_low = rep(NA,n_weight)
    auc = rep(NA,n_weight)
    auc_high = rep(NA,n_weight)
    x = vector("list",n_weight)
    y = vector("list",n_weight)
    for (i_weight in 1:n_weight) {
        pred = c(as.numeric(pval_mat[[1]][,paste(type,weight[i_weight],sep="_")]),
        as.numeric(pval_mat[[2]][,paste(type,weight[i_weight],sep="_")]))
        res = roc(true~pred,ci=T,direction=">")
        # direction ">" if the predictor (p-value) is higher for controls than cases
        # That is, here we are explicitly stating that positive controls should have a lower p-value
        # we compute the 95% CI (DeLong)
        auc_low[i_weight] = round(ci.auc(res)[1],digits=2)
        auc[i_weight] = round(ci.auc(res)[2],digits=2)
        auc_high[i_weight] = round(ci.auc(res)[3],digits=2)
        y[[i_weight]] = res$sensitivities # Sensitivity (True Positive Rate)
        x[[i_weight]] = 1-res$specificities # 1 - Specificity (False Positive Rate)
    }
    title = paste0(proj[i_proj],"(",type,"): ",collection[1]," vs ",collection[2])
    plot(c(0,1),c(0,1),type="n",xlab="1 - Specificity (False Positive Rate)",
    ylab="Sensitivity (True Positive Rate)",cex.lab=1,main=title)
    col = 1:n_weight
    for (i_weight in 1:n_weight) {
        points(x[[i_weight]],y[[i_weight]],
        cex=1,col=col[i_weight],lty=1,type="l",lwd=2)
    }
    points(c(0,1),c(0,1),type="l",lty=2) # ref line
    legend = rep(NA,n_weight)
    for (i_weight in 1:n_weight) {
        legend[i_weight] = paste0(weight[i_weight],": AUC=",auc[i_weight]," (95% CI:",auc_low[i_weight],"-",auc_high[i_weight],")")
    }
    legend("bottomright",legend,col=col,lty=1,lwd=2,cex=0.85)
    dev.off()
    
    outfile = file.path(output_dir,paste0(proj[i_proj],"_",type,"_AUC_",paste0(collection,collapse="_"),".txt"))
    res = cbind(weight,auc_low,auc,auc_high)
    output = rbind(colnames(res),res)
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
