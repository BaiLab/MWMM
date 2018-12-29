rm(list=ls())
options(stringsAsFactors = F)
ptm <- proc.time()
args<-commandArgs(TRUE)
library("magrittr")
# source("https://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")
# source("https://bioconductor.org/biocLite.R")
# biocLite("DOSE")
# biocLite("GOSemSim")

# -----------------------------------------------------------------------------
source(file.path(getwd(), "function_scripts_directory", "calculate_enrichment_from_a_list_of_cluster.R"))
source(file.path(getwd(), "function_scripts_directory", "GO_similarity_functions.R"))
#------------------source functions below---------------------------------------
conversion_table_path.TCGA <- file.path(getwd(), "function_scripts_directory", "mRNA_names_and_Entrez_IDs_dictionary.TCGA.csv")

conversion_table_path.miRMAP_bicluster <- file.path(getwd(), "function_scripts_directory", "mRNA_names_and_Entrez_IDs_dictionary.miRMAP_bicluster.csv")


TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])

cancer_types <- rev(cancer_types)

for (cancer_type in cancer_types) {
    
    output_directory <- file.path(getwd(),
                                  "output_data_directory",
                                  "binary_object_files",
                                  paste("enrichment_analysis", cancer_type, sep = "."))
    dir.create(output_directory, showWarnings = F, recursive = T)
    
    
    # get algorithm names, ie, communities_rounds.
    path_to_communities_input_directory <- file.path(getwd(),
                                                     "output_data_directory",
                                                     "binary_object_files",
                                                     paste("communities_files", cancer_type, sep = "."))
    
    
    # get the names of different algorithms.
    communities_rounds <- dir(path_to_communities_input_directory) %>%
        strsplit(., "\\.") %>%
        sapply(., function(x) x[2][[1]])
    
    
    input_file_names <- paste("the_communities", communities_rounds, "rds", sep = ".")
    path_to_input_files <- file.path(getwd(),
                                     "output_data_directory",
                                     "binary_object_files",
                                     paste("communities_files", cancer_type, sep = "."),
                                     input_file_names)
    
    
    a_list_of_the_communities <- lapply(path_to_input_files, readRDS)
    names(a_list_of_the_communities) <- communities_rounds
    
    
    
    # -------------convert clusters in communities object to a list of list of clusters
    a_list_of_cluster_list <- vector(mode = "list", length = length(communities_rounds))
    
    for (i in seq_along(communities_rounds)) {
        path_to_input_file <- path_to_input_files[i]
        current_communities <- readRDS(path_to_input_file)
        
        if (communities_rounds[i] == "miRMAP_bicluster") {
            current_cluster_list <- convert_a_communities_to_a_list_of_mRNA_only_clusters(current_communities, conversion_table_path.miRMAP_bicluster)
        } else {
            current_cluster_list <- convert_a_communities_to_a_list_of_mRNA_only_clusters(current_communities, conversion_table_path.TCGA)
        }
        
        # filter out the clusters that only have miRNAs.
        current_cluster_list <- current_cluster_list[!sapply(current_cluster_list, is.list)]
        
        a_list_of_cluster_list[[i]] <- current_cluster_list
        # path_to_output_file <- file.path(getwd(), "output_data_directory", "plain_text_files", "GO_similarity_distance.BRCA", output_file_names[i])
        # write.csv(GO_similarity_df, file = path_to_output_file, row.names = F, quote = F)
        print(i)
    }
    names(a_list_of_cluster_list) <- communities_rounds
    
    
    # all_indicator_list <- foreach(algorithm_name = names(a_list_of_cluster_list)) %dopar% {
    # the_cluster_list <- a_list_of_cluster_list$blossom_02
    
    for (algorithm_name in names(a_list_of_cluster_list)) {
        this_cluster_list <- a_list_of_cluster_list[[algorithm_name]]
        
        
        current_enriched_terms_list <- calculate_enrichment_from_a_cluster_list(this_cluster_list)
        output_name <- paste("biological_term_list", algorithm_name, "rds", sep = ".")
        output_object_path <- file.path(output_directory, output_name)
        saveRDS(current_enriched_terms_list, file = output_object_path)
        
        cat(cancer_type, algorithm_name, "is done", sep = " ", "\n")
    }
}




proc.time() - ptm



# a_list_of_cluster_list$blossom_02$`38`

# --------------run the biological term analysis and return the result to a list.
# all_indicator_list <- list(mode="list", length = length(a_list_of_cluster_list))

# calculate_enrichment_from_a_list_of_cluster(this_cluster_list)
# # gene-disease associations
# ck_trial <- try(compareCluster(geneCluster = the_cluster_list,
#                                fun = "enrichDGN",
#                                pvalueCutoff  = 0.01,
#                                pAdjustMethod = "BH",
#                                minGSSize = 10,
#                                maxGSSize = 500,
#                                qvalueCutoff = 0.05))
# 
# if (is.null(attr(ck_trial, "condition")$message)) {
#     ck <- ck_trial
#     print("ck is done, because there is no error")
# } else {
#     ck <- NULL
#     print( "ck is not done, because no enrichment of enrichDGN")
# }
# 
# # Disease-Ontology
# # do <- compareCluster(geneClusters = the_cluster_list, fun = "enrichDO")
# 
# do_trial <- try(compareCluster(geneClusters = the_cluster_list,
#                                fun = "enrichDO",
#                                ont = "DO",
#                                pvalueCutoff  = 0.01,
#                                pAdjustMethod = "BH",
#                                minGSSize = 10,
#                                maxGSSize = 500,
#                                qvalueCutoff = 0.05,
#                                readable = FALSE))
# 
# if (is.null(attr(do_trial, "condition")$message)) {
#     do <- do_trial
#     print("do is done, because there is no error")
# } else {
#     do <- NULL
#     print( "do is not done, because no enrichment of enrichDO")
# }
# 
# 
# # Network-of-Cancer-Gene
# ncg_trial <- try(compareCluster(geneClusters = the_cluster_list,
#                                 fun = "enrichNCG",
#                                 pvalueCutoff  = 0.01,
#                                 pAdjustMethod = "BH",
#                                 minGSSize = 10,
#                                 maxGSSize = 500,
#                                 qvalueCutoff = 0.05))
# 
# if (is.null(attr(ncg_trial, "condition")$message)) {
#     # no error message pops up.
#     ncg <- ncg_trial
#     print("ncg is done, because there is no error.")
# } else {
#     ncg <- NULL
#     print( "ncg is not done, because no enrichment of enrichDO")
# }
# 
# # kegg-over-representation
# kk_trial <- try(compareCluster(geneClusters = the_cluster_list,
#                                fun = "enrichKEGG",
#                                pvalueCutoff  = 0.01,
#                                pAdjustMethod = "BH",
#                                minGSSize = 10,
#                                maxGSSize = 500,
#                                qvalueCutoff = 0.05))
# 
# if (is.null(attr(kk_trial, "condition")$message)) {
#     # no error message pops up.
#     kk <- kk_trial
#     print("kk is done, because there is no error.")
# } else {
#     kk <- NULL
#     print( "kk is not done, because no enrichment of enrichKEGG")
# }
# # kk@geneClusters[[38]]
# # kk_df <- kk@compareClusterResult
# 
# 
# # GO-over-representation
# ontology_words <- c("BP", "CC", "MF")
# # enriched_GO_list <- list(mode="list", length = length(ontology_words))
# 
# enriched_GO_list <- foreach(j = seq_along(ontology_words)) %dopar% {
#     # for (j in seq_along(ontology_words)) {
#     ego_trial <- try(compareCluster(geneClusters = the_cluster_list,
#                                     fun = "enrichGO",
#                                     OrgDb = org.Hs.eg.db,
#                                     ont = ontology_words[j],
#                                     pvalueCutoff  = 0.01,
#                                     pAdjustMethod = "BH",
#                                     minGSSize = 10,
#                                     maxGSSize = 500,
#                                     qvalueCutoff = 0.05,
#                                     readable= FALSE))
#     
#     if (is.null(attr(ego_trial, "condition")$message)) {
#         # no error message pops up.
#         ego <- ego_trial
#         print("ego is done, because there is no error.")
#     } else {
#         ego <- NULL
#         print( "ego is not done, because no enrichment of enrichGO")
#     }
#     # enriched_GO_list[[j]] <- ego
#     ego
# }
# 
# names(enriched_GO_list) <- ontology_words
# ego_BP <- enriched_GO_list[["BP"]]
# ego_CC <- enriched_GO_list[["CC"]]
# ego_MF <- enriched_GO_list[["MF"]]
# 
# print("GO is done")
# 
# # ------------------compile the biological terms into a list.------
# current_indicator_list <- list(ck, do, ncg, kk, ego_BP, ego_CC, ego_MF)
# names(current_indicator_list) <- c("ck", "do", "ncg", "kk", "ego_BP", "ego_CC", "ego_MF")
# 
# # -----------------save and return results.-------------------------
# output_name <- paste("biological_term_list", algorithm_name, "rds", sep = ".")
# output_object_path <- file.path(output_directory, output_name)
# saveRDS(current_indicator_list, file = output_object_path)
# 
# all_indicator_list[[algorithm_name]] <- current_indicator_list
# 
# cat(algorithm_name, "is done", sep = " ", "\n")
# # current_indicator_list
# }
# 
# 
# names(all_indicator_list) <- communities_rounds
# 
# # path_to_output_file <- file.path(getwd(), "output_data_directory", "binary_object_files", "enrichment_analysis.BRCA", "all_indicator_list.rds")
# # saveRDS(all_indicator_list, path_to_output_file)
# 
# 
# stopCluster(current_cluster)

# x <- readRDS(file.path(getwd(), "output_data_directory", "binary_object_files", "enrichment_analysis.BRCA", "biological_term_list.leading_eigen.rds"))
# 
# x_ck <- x[["ck"]]
# y <- x_ck@compareClusterResult
# 
# x_go_cc <- x$ego_CC@compareClusterResult
# x_go_bp <- x$ego_BP@compareClusterResult
# x_go_mf <- x$ego_MF@compareClusterResult


# x <- org.Hs.egGO
# hsEG <- mappedkeys(x)
# set.seed <- 123
# sim_db <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
# clusters <- list(a=sample(hsEG, 20), b=sample(hsEG, 20), c=sample(hsEG, 20))
# y <- mclusterSim(clusters, semData = sim_db, measure="Wang", combine="BMA")
# combineScores(y, combine = "avg")

#-------------------prepare input files--------------
# communities_rounds <- c("label_propagation",
#                         "walktrap_algorithm",
#                         "fast_greedy",
#                         "leading_eigen",
#                         "edge_betweenness",
#                         "markov_matrix",
#                         "louvain_algorithm",
#                         "miRMAP_bicluster",
#                         "hungarian_algorithm",
#                         "blossom_01",
#                         "blossom_02",
#                         "blossom_03",
#                         "blossom_04",
#                         "blossom_05",
#                         "blossom_06",
#                         "blossom_07",
#                         "blossom_08")

# command line take the algorithm name index.
# communities_rounds <- communities_rounds[as.numeric(args[1])]
# communities_rounds <- rev(communities_rounds)

# cat("current communities is", communities_rounds, sep=" ", "\n")
