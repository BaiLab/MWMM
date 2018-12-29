rm(list=ls())
options(stringsAsFactors = F)
# options(device = "pdf")
# getOption("device")

ptm <- proc.time()
args<-commandArgs(TRUE)

library("magrittr")
library("clusterProfiler")
library("ggplot2")

# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])


output_graph_directory <- file.path(getwd(), "output_data_directory", "graph_files", "cluster_enrichment_graphs")
dir.create(output_graph_directory, showWarnings = F, recursive = T)

for (cancer_type in cancer_types) {
    # cancer_type <- "LUAD"
    
    # get algorithm names, ie, communities_rounds.
    path_to_communities_input_directory <- file.path(getwd(),
                                                     "output_data_directory",
                                                     "binary_object_files",
                                                     paste("communities_files", cancer_type, sep = "."))
    
    
    # get the names of different algorithms.
    communities_rounds <- dir(path_to_communities_input_directory) %>%
        strsplit(., "\\.") %>%
        sapply(., function(x) x[2][[1]])
    
    our_communities_rounds <- c("hungarian_algorithm", communities_rounds[grepl("blossom", communities_rounds)])
    other_communities_rounds <- c("fast_greedy",
                                  "leading_eigen",
                                  "edge_betweenness",
                                  "label_propagation",
                                  "louvain_algorithm",
                                  "walktrap_algorithm")
    
    communities_rounds <- c(other_communities_rounds, our_communities_rounds)
    
    
    input_file_names <- paste("biological_term_list", communities_rounds, "rds", sep = ".")
    path_to_input_files <- file.path(getwd(),
                                     "output_data_directory",
                                     "binary_object_files",
                                     paste("enrichment_analysis", cancer_type, sep = "."),
                                     input_file_names)
    
    a_list_of_enrichment_list <- lapply(path_to_input_files, readRDS)
    names(a_list_of_enrichment_list) <- communities_rounds
    
    
    for (i in names(a_list_of_enrichment_list)) {
        # i <- "hungarian_algorithm"
        
        an_enrichment_result_list <- a_list_of_enrichment_list[[i]]
        
        #  an enrichment result of a clustering algorithm contains
        # "ck", "do", "ncg", "kk", "ego_BP", "ego_CC", "ego_MF"
        for (j in names(an_enrichment_result_list)) {
            
            # j <- "kk"
            a_term_object <- an_enrichment_result_list[[j]]
            if (is.null(a_term_object)) {
                next
            } else {
                output_file_name <- paste("enrichment", cancer_type, i, j, "pdf", sep = ".")
                path_to_graph_output <- file.path(output_graph_directory, output_file_name)
                
                # The dotplot is calling ggplot2. If you understand Chinese, please
                # refer to https://guangchuangyu.github.io/cn/2017/07/clusterprofiler-dotplot/
                
                if (j == "kk") {
                    y_axis_label <- "Enriched pathways"
                } else {
                    y_axis_label <- "Enriched terms"
                }
                
                p <- dotplot(a_term_object)
                q <- p + labs(y=y_axis_label, x = "Clusters")
                ggsave(path_to_graph_output, q, width = 170, height = 170, units = "mm", scale = 2.5)
                # ggsave(path_to_graph_output, q, width = 20, height = 20, units = "in")
            }
        }
        
    }
    # break
}

proc.time() - ptm

# 
# #-------------------prepare input files--------------
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
# 
# # command line take the algorithm name index.
# 
# 
# cat("current communities is", communities_rounds, sep=" ", "\n")
# 
# 
# input_file_names <- paste("biological_term_list", communities_rounds, "rds", sep = ".")
# path_to_input_files <- file.path(getwd(), "output_data_directory", "binary_object_files", "enrichment_analysis.BRCA", input_file_names)
# 
# a_list_of_enrichment_list <- lapply(path_to_input_files, readRDS)
# names(a_list_of_enrichment_list) <- communities_rounds
# 
# 
# # extract the biological terms from the list and calculate the average.
# accum_list <- vector(mode = "list", length = length(communities_rounds))
# for (i in names(a_list_of_enrichment_list)) {
#     
#     # i <- "hungarian_algorithm"
#         
#     an_enrichment_result_list <- a_list_of_enrichment_list[[i]]
#     
#     current_terms_vector <- numeric(length(an_enrichment_result_list))
#     names(current_terms_vector) <- names(an_enrichment_result_list)
#     
#     #  an enrichment result of a clustering algorithm contains
#     # "ck", "do", "ncg", "kk", "ego_BP", "ego_CC", "ego_MF"
#     for (j in names(an_enrichment_result_list)) {
#         
#         # j <- "kk"
#         a_term_object <- an_enrichment_result_list[[j]]
#         if (is.null(a_term_object)) {
#             mean_terms_clusters <- NA
#         } else {
#             a_term_df <- a_term_object@compareClusterResult
#             
#             out_path <- file.path(getwd(),
#                                   "output_data_directory",
#                                   "plain_text_files",
#                                   "biological_terms.BRCA")
#             
#             output_name <- paste(i, j, "tsv", sep = ".")
#             write.table(a_term_df, file = file.path(out_path, output_name),
#                       quote = F, row.names = F, sep = "\t")
#             
#             cluster_number <- length(a_term_object@geneClusters)
#             
#             term_frequency_cluster <- table(a_term_df$Cluster)
#             
#             cluster_length_vector <- sapply(a_term_object@geneClusters, length)
#             
#             existing_term_cluster <- ifelse(term_frequency_cluster != 0, 1, 0)
#             
#             # Three ways of calculating enrichment scores of different algorithm.
#             
#             mean_terms_clusters <- mean(term_frequency_cluster / cluster_length_vector)
#             mean_terms_clusters <- sum(term_frequency_cluster != 0)/cluster_number
#             mean_terms_clusters <- mean(existing_term_cluster / cluster_length_vector)
#             
#         }
#         current_terms_vector[j] <- mean_terms_clusters
#         print(mean_terms_clusters)
#         print("_______________________________")
#     }
#     accum_list[[i]] <- current_terms_vector
# }
# 
# algorithm_name <- "blossom_01"
# biological_term <- "kk"
# x_term_object <- a_list_of_enrichment_list[[algorithm_name]][[biological_term]]
# x_term_df <- x_term_object@compareClusterResult
# 
# x_term_object@geneClusters[[52]]
# 
# 
# dotplot(x_term_object)
# 
# # cluster validation
# 
# input_file_names <- paste("the_communities", communities_rounds, "rds", sep = ".")
# path_to_input_files <- file.path(getwd(), "output_data_directory", "binary_object_files", "communities_files.BRCA", input_file_names)
# 
# 
# a_list_of_the_communities <- lapply(path_to_input_files, readRDS)
# names(a_list_of_the_communities) <- communities_rounds
# x <- a_list_of_the_communities[["blossom_01"]]
# 
# communities(x)[52]
# pdf("plots.pdf")
# dev.off()
# 
# 
# file_type <- "pdf"
# output_graph_file_name <- paste("KEGG", algorithm_name, file_type, sep = ".")
# 
# ggsave(output_graph_file_name, plot = last_plot(),
#        width = 170, height = 170, units = "mm", device = file_type,
#        path = file.path(getwd(), "output_data_directory", "graph_files"))

# final_df <- do.call(rbind, accum_list)
# final_df <- round(final_df, digits = 3)
# 
# new <- a_list_of_enrichment_list$hungarian_algorithm
# cluster_x <- new$ego_BP@geneClusters[[186]]
# length(cluster_x)




