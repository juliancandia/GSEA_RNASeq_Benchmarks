library(pROC)

rm(list=ls())

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

#type = "signed"
type = "unsigned"

method = c("Bonf","BH")
method_label = c("Bonf","B-H")
n_method = length(method)

#collection = c("pwy.presel","neg.ctrl.unif") # (TRUE,NEG) collections to compare
collection = c("pos.ctrl","neg.ctrl.unif") # (TRUE,NEG) collections to compare
n_coll = 2

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

for (i_proj in 1:n_proj) {
    output_dir = file.path(PROJECT_DIR,proj[i_proj])
    pval_mat = vector("list",n_coll)
    for (i_coll in 1:n_coll) {
        for (i_method in 1:n_method) {
            infile = file.path(output_dir,paste0(proj[i_proj],"_",collection[i_coll],"_ORA_",method[i_method],".txt"))
            data = as.matrix(read.table(infile,header=T,stringsAsFactors=F,sep="\t"))
            if (type=="signed") {
                pval_mat[[i_coll]] = cbind(pval_mat[[i_coll]],as.numeric(data[,"pval_merged"]))
            } else if (type=="unsigned") {
                pval_mat[[i_coll]] = cbind(pval_mat[[i_coll]],as.numeric(data[,paste0("pval_abs_",method[i_method])]))
            }
        }
    }
    
    # we generate ROC plots / AUC results
    outfile = file.path(output_dir,paste0(proj[i_proj],"_ORA_",type,"_ROC_",paste0(collection,collapse="_"),".pdf"))
    pdf(outfile,width=7,height=5)
    true = c(rep(1,nrow(pval_mat[[1]])),rep(0,nrow(pval_mat[[2]])))
    auc_low = rep(NA,n_method)
    auc = rep(NA,n_method)
    auc_high = rep(NA,n_method)
    x = vector("list",n_method)
    y = vector("list",n_method)
    for (i_method in 1:n_method) {
        pred = c(as.numeric(pval_mat[[1]][,i_method]),as.numeric(pval_mat[[2]][,i_method]))
        res = roc(true~pred,ci=T,direction=">")
        # direction ">" if the predictor (p-value) is higher for controls than cases
        # That is, here we are explicitly stating that positive controls should have a lower p-value
        # we compute the 95% CI (DeLong)
        auc_low[i_method] = round(ci.auc(res)[1],digits=2)
        auc[i_method] = round(ci.auc(res)[2],digits=2)
        auc_high[i_method] = round(ci.auc(res)[3],digits=2)
        y[[i_method]] = res$sensitivities # Sensitivity (True Positive Rate)
        x[[i_method]] = 1-res$specificities # 1 - Specificity (False Positive Rate)
    }
    title = paste0(proj[i_proj]," (",type," ORA): ",collection[1]," vs ",collection[2])
    plot(c(0,1),c(0,1),type="n",xlab="1 - Specificity (False Positive Rate)",
    ylab="Sensitivity (True Positive Rate)",cex.lab=1,main=title)
    col = 1:n_method
    for (i_method in 1:n_method) {
        points(x[[i_method]],y[[i_method]],
        cex=1,col=col[i_method],lty=1,type="l",lwd=2)
    }
    points(c(0,1),c(0,1),type="l",lty=2) # ref line
    legend = rep(NA,n_method)
    for (i_method in 1:n_method) {
        legend[i_method] = paste0(method_label[i_method],": AUC=",auc[i_method]," (95% CI:",auc_low[i_method],"-",auc_high[i_method],")")
    }
    legend("bottomright",legend,col=col,lty=1,lwd=2,cex=0.85)
    dev.off()
    
    outfile = file.path(output_dir,paste0(proj[i_proj],"_ORA_",type,"_AUC_",paste0(collection,collapse="_"),".txt"))
    res = cbind(method,auc_low,auc,auc_high)
    output = rbind(colnames(res),res)
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
