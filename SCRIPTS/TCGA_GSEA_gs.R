rm(list=ls())

# programmatic approach to GSEA following Reimand et al, Nat Prot 2019

PROJECT_DIR = "~/NIA/PATHWAYS/BENCHMARK/GitHub/RESULTS" # replace with your local path

collection = "pwy.presel" # pathways preselected as positive controls
# collection = "neg.ctrl.unif" # randomized negative control pathways with uniform distribution of pathway sizes between size_min=15 and size_max=500

weight_label = c("cl",paste0("p",c("1","1.5","2")))
weight = c("classic",paste0("weighted",c("","_p1.5","_p2")))
n_weight = length(weight)

infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list.txt")
#infile = file.path(PROJECT_DIR,"TCGA","TCGA_proj_list_filt.txt")
proj = as.matrix(read.table(infile,header=F,sep="\t",stringsAsFactors=F))[,1]
n_proj = length(proj)

gsea_script = file.path("~/NIA/PATHWAYS","GSEA","GSEA_4.3.2","gsea-cli.sh") # replace with your local path
for (i_proj in 1:n_proj) {
    # GSEA Preranked
    # User Guide: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Other_Ways_to
    output_dir = file.path(PROJECT_DIR,proj[i_proj])
    log_dir = file.path(output_dir,"log")
    dir.create(log_dir,showWarnings=F)
    if (collection=="pwy.presel") {
        gmt_file = file.path(output_dir,paste0(proj[i_proj],"_pwy_presel.gmt"))
    } else if (collection=="neg.ctrl.unif") {
        gmt_file = file.path(output_dir,paste0(proj[i_proj],"_pwy_neg.ctrl.unif.gmt"))
    }
    rnk_file = file.path(output_dir,paste0(proj[i_proj],"_ranks.rnk"))
    n_perm = 1000
    seed = "timestamp"
    n_plot = 2
    size_max = 500
    size_min = 15
    for (i_weight in 1:n_weight) {
        run_label = paste0(proj[i_proj],"_",collection,"_gs_",weight_label[i_weight])
        command <- paste(gsea_script,"GSEAPreranked -gmx",gmt_file,"-collapse No_Collapse -mode Abs_max_of_probes -norm meandiv -nperm",n_perm,"-rnd_seed",seed,
        "-rnk",rnk_file,
        "-scoring_scheme",weight[i_weight],"-rpt_label",run_label,"-create_svgs false -include_only_symbols true -make_sets true -plot_top_x",n_plot,"-set_max",size_max,"-set_min",size_min," -zip_report false -out",output_dir,">", file.path(log_dir,paste0("log_",run_label,".txt"))
        )
        system(command)
    }
}
