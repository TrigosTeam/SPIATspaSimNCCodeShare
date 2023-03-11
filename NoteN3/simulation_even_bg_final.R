# This script is for generating simulated images.

load("Objects/independent_sim_args.Rda") # This object was generated from "diabetes_extract_args.R"
load("Objects/even_bg_paired.Rda") # This object was generated from "simulate_even_bg_cells.R"

args <- args_sim_list$args

args <-merge(args, args_sim_list$infil_idx_t)
args <-merge(args, args_sim_list$stromal_idx_t)
args <- merge(args, args_sim_list$dist_idx_t)

# calculate dependent parameters ####
# Image size should be extracted from the simulated background images
for (i in 1:length(even_bg_paired)){
    image <- even_bg_paired[[i]]$image
    args[i, "X_max"] <- max(SpatialExperiment::spatialCoords(image)[, "Cell.X.Position"]) 
    args[i, "Y_max"] <- max(SpatialExperiment::spatialCoords(image)[, "Cell.Y.Position"]) 
}
# function calculate_cluster_size is updated! Using the density to calculate the cluster size
calculate_cluster_size <- function(args){
    # if there are 4, 5 or 6 clusters, treat it as 4 clusters
    args[, "n_clusters2"] <- as.numeric(args[, "n_clusters"])
    args[which(args[, "n_clusters2"] > 4), "n_clusters2"] <- 4
    # calculate the number of cluster cells in each clusters
    args[, "n_cluster_cells2"] <- as.numeric(args[, "n_cluster_cells"])/as.numeric(args[, "n_clusters2"])
    # ## This should use average minimum dist, WAS CORRECT
    # # The previous simulation using average min dist was wrong, should have used min min dist
    # args[, "cluster_size2"] <-
    #     ceiling(2* sqrt(args[, "n_cluster_cells2"]) * as.numeric(args[, "avg_min_dist"])/sqrt(pi))
    # use the density to calculate the radius
    args[, "cluster_size"] <- sqrt(args[, "n_cluster_cells2"]/(0.009401827*pi))
    return(args)
}

new <- calculate_cluster_size(args)


calculate_image_size <- function(args){
    # # This is derived from the ANN calculation
    # min_dist <- as.numeric(args[, "avg_min_dist"])
    n_cells <- as.numeric(args[, "n_cells"])
    # calculate the area 
    # A <- 4 * min_dist**2 * n_cells
    A<- n_cells/0.009401827
    args[, "X_max"] <- round(sqrt(A/as.numeric(args[, "window_ratio"])), 0)
    args[, "Y_max"] <- round(args[, "X_max"] * as.numeric(args[, "window_ratio"]), 0)
    return(args)
}

new <- calculate_image_size(new)

calculate_centre_locations <- function(x){
    for (i in seq_len(nrow(x))){
        if (x[i, "n_clusters2"] == 1) {
            x[i, "cluster_1_x"] <- as.numeric(x[i, "X_max"])/2
            x[i, "cluster_1_y"] <- as.numeric(x[i, "Y_max"])/2
        }
        else if(x[i, "n_clusters2"] == 2) {
            x[i, "cluster_1_x"] <- as.numeric(x[i, "X_max"])/3
            x[i, "cluster_1_y"] <- as.numeric(x[i, "Y_max"])/3 *2
            
            x[i, "cluster_2_x"] <- as.numeric(x[i, "X_max"])/3 * 2
            x[i, "cluster_2_y"] <- as.numeric(x[i, "Y_max"])/3
        }
        else if(x[i, "n_clusters2"] == 3) {
            x[i, "cluster_1_x"] <- as.numeric(x[i, "X_max"])/2
            x[i, "cluster_1_y"] <- as.numeric(x[i, "Y_max"])/3 *2
            
            x[i, "cluster_2_x"] <- as.numeric(x[i, "X_max"])/3
            x[i, "cluster_2_y"] <- as.numeric(x[i, "Y_max"])/3
            
            x[i, "cluster_3_x"] <- as.numeric(x[i, "X_max"])/3 * 2
            x[i, "cluster_3_y"] <- as.numeric(x[i, "Y_max"])/3
        }
        else{
            x[i, "cluster_1_x"] <- as.numeric(x[i, "X_max"])/4
            x[i, "cluster_1_y"] <- as.numeric(x[i, "Y_max"])/4
            
            x[i, "cluster_2_x"] <- as.numeric(x[i, "X_max"])/3 * 2
            x[i, "cluster_2_y"] <- as.numeric(x[i, "Y_max"])/3
            
            x[i, "cluster_3_x"] <- as.numeric(x[i, "X_max"])/4 * 3
            x[i, "cluster_3_y"] <- as.numeric(x[i, "Y_max"])/4 * 3
            
            x[i, "cluster_4_x"] <- as.numeric(x[i, "X_max"])/3
            x[i, "cluster_4_y"] <- as.numeric(x[i, "Y_max"])/3 * 2
        }
    }
    return(x)
}

new <- calculate_centre_locations(new)

library(spaSim)
library(SPIAT)

even_sim_paired <- list()

infiltration_types <- c("nonbeta", "immune", "endothelial", "others")
for (i in 1:845){
    args <- new[i, ]

    bg <- even_bg_paired[[i]]$image
    
    mixing_bg <- simulate_mixing(
        bg_sample = bg, idents = c("beta", infiltration_types),
        props = c(
            as.numeric(args$p_stromal_beta),as.numeric(args$p_stromal_nonbeta),
            as.numeric(args$p_stromal_immune),as.numeric(args$p_stromal_endothelial),
            as.numeric(args$p_stromal_others)
        ),plot_image = FALSE)
    
    if (as.numeric(args$n_clusters2) == 1){
        image <- simulate_clusters(bg_sample = mixing_bg, n_clusters = 1, cluster_properties = list(
            C1 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_1_x), y = as.numeric(args$cluster_1_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            )
        ), plot_image = FALSE)
    }else if (as.numeric(args$n_clusters2) == 2){
        image <- simulate_clusters(bg_sample = mixing_bg, n_clusters = 2, cluster_properties = list(
            C1 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_1_x), y = as.numeric(args$cluster_1_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            ),
            C2 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_2_x), y = as.numeric(args$cluster_2_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            )), plot_image = FALSE)
    }
    else if (as.numeric(args$n_clusters2) == 3){
        image <- simulate_clusters(bg_sample = mixing_bg, n_clusters = 3, cluster_properties = list(
            C1 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_1_x), y = as.numeric(args$cluster_1_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            ),
            C2 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_2_x), y = as.numeric(args$cluster_2_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            ),
            C3 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_3_x), y = as.numeric(args$cluster_3_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            )), plot_image = FALSE)
    }else{
        image <- simulate_clusters(bg_sample = mixing_bg, n_clusters = 4, cluster_properties = list(
            C1 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_1_x), y = as.numeric(args$cluster_1_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            ),
            C2 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_2_x), y = as.numeric(args$cluster_2_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            ),
            C3 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_3_x), y = as.numeric(args$cluster_3_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            ),
            C4 = list(
                name_of_cluster_cell = "beta",
                size = as.numeric(args$cluster_size),
                shape = "Circle",
                centre_loc = data.frame(x = as.numeric(args$cluster_4_x), y = as.numeric(args$cluster_4_y)),
                infiltration_types = infiltration_types,
                infiltration_proportions = c(
                    as.numeric(args$p_infil_nonbeta), as.numeric(args$p_infil_immune), 
                    as.numeric(args$p_infil_endothelial), as.numeric(args$p_infil_others)
                )
            )), plot_image = FALSE)
    }
    # }
    image_spe <- format_colData_to_spe(image)
    attr(image_spe, "name") <- paste0("Simulation_", i)
    g <- plot_cell_categories(image_spe, categories_of_interest = c("beta", infiltration_types),
                              colour_vector = c("orange", "darkgreen", "red", "blue", "gray"), 
                              feature_colname = "Cell.Type")
    print(g)
    even_sim_paired[[i]] <- image_spe
    save(even_sim_paired, file = "Objects/even_sim_paired.Rdata")
}


even_sim_paired_final <- list()
n_cluster_cells_sim <- c()
for (i in 1:845){
    spe <- even_sim_paired[[i]]
    spe <- define_celltypes(spe, categories = c("beta", "nonbeta", "immune", 
                                                "endothelial", "others"), 
                            category_colname = "Cell.Type",
                            names = c("islet", "islet", "immune", "endothelial",
                                      "others"), new_colname = "Cell.Type2")
    spe_border <- identify_bordering_cells(spe, reference_cell = "islet", 
                                           feature_colname = "Cell.Type2", 
                                           ahull_alpha = 40)
    counts <- table(spe_border@colData@listData[["Region"]])
    n_cluster_cells <- counts["Border"] + counts["Inside"] 
    n_cluster_cells_sim <- c(n_cluster_cells_sim, n_cluster_cells)
    colData(spe_border)$Phenotype <- NULL
    even_sim_paired_final[[i]] <- spe_border
    print(i)
    save(even_sim_paired_final, file = "Objects/even_sim_paired_final.Rda")
}
