load("Objects/diabetes_real_spe.Rda")
load("Objects/even_sim_paired_final.Rda")

final_metrics <- data.frame(ID = 1:845)

# get the number of total cells per image and window size ####
n_cells_real <- c()
n_cells_sim <- c()
width_real <- c()
width_sim <- c()
height_real <- c()
height_sim <- c()

for (i in 1:845){
    image_real <- Real[[i]]
    image_sim <- even_sim_paired_final[[i]]
    n_cells_real <- c(n_cells_real, ncol(image_real))
    n_cells_sim <- c(n_cells_sim, ncol(image_sim))
    
    width_real <- c(width_real, max(spatialCoords(image_real)[, 1]))
    width_sim <- c(width_sim, max(spatialCoords(image_sim)[, 1]))
    height_real <- c(height_real, max(spatialCoords(image_real)[, 2]))
    height_sim <- c(height_sim, max(spatialCoords(image_sim)[, 2]))
}

final_metrics <- cbind(final_metrics, n_cells_real, n_cells_sim, 
                       width_real, width_sim, height_real, 
                       height_sim)

# get the clustering metrics of each image ####
n_cluster_cells_real <- c()
n_cluster_cells_sim <- c()
R_BC_real <- c()
R_BC_sim <- c()
for (i in 1:845){
    image_real <- Real[[i]]
    image_sim <- even_sim_paired_final[[i]]
    n_cluster_cells_real <- c(n_cluster_cells_real, length(which(image_real$Region %in% c("Inside", "Border"))))
    n_cluster_cells_sim <- c(n_cluster_cells_sim, length(which(image_sim$Region %in% c("Inside", "Border"))))
    R_BC_real <- c(R_BC_real, R_BC(image_real, cell_type_of_interest = "islet", feature_colname = "Cell.Type4"))
    R_BC_sim <- c(R_BC_sim, R_BC(image_sim, cell_type_of_interest = "islet", feature_colname = "Cell.Type2"))
}

final_metrics <- cbind(final_metrics, n_cluster_cells_real, n_cluster_cells_sim,
                       R_BC_real, R_BC_sim)
save(final_metrics, file = "Objects/final_metrics.Rda")

# get the colocalisation metrics between islet and endothelial cells of each image ####
ms_real <- c()
ms_sim <- c()
nms_real <- c()
nms_sim <- c()
auc_real <- c()
auc_sim <- c()
avg_p_dist_real <- c()
avg_p_dist_sim <- c()
avg_m_dist1_real <- c()
avg_m_dist2_real <- c()
avg_m_dist1_sim <- c()
avg_m_dist2_sim <- c()
cin_real <- c()
cin_sim <- c()

for (i in 1:845){
    image_real <- Real[[i]]
    image_sim <- even_sim_paired_final[[i]]
    # ms & nms
    radius_real <- 3 * average_minimum_distance(image_real)
    ms1 <- mixing_score_summary(image_real, reference_celltype = "islet", 
                                target_celltype = "endothelial", 
                                feature_colname = "Cell.Type4",
                                radius = radius_real)
    radius_sim <- 3 * average_minimum_distance(image_sim)
    ms2 <- mixing_score_summary(image_sim, reference_celltype = "islet", 
                                target_celltype = "endothelial", 
                                feature_colname = "Cell.Type2",
                                radius = radius_sim)
    ms_real <- c(ms_real, ms1[, "Mixing_score"])
    nms_real <- c(nms_real, ms1[, "Normalised_mixing_score"])
    ms_sim <- c(ms_sim, ms2[, "Mixing_score"])
    nms_sim <- c(nms_sim, ms2[, "Normalised_mixing_score"])
    
    # AUC 
    crossK <- try(calculate_cross_functions(image_real, cell_types_of_interest = c("islet", "endothelial"),
                                               feature_colname = "Cell.Type4"))
    if("try-error" %in% class(crossK) ){
        auc  <- NA
    } else{
        auc <- AUC_of_cross_function(crossK) 
    }
    
    auc_real <- c(auc_real, auc)
    
    crossK <- try(calculate_cross_functions(image_sim, cell_types_of_interest = c("islet", "endothelial"),
                                            feature_colname = "Cell.Type2"))
    if("try-error" %in% class(crossK) ){
        auc  <- NA
    } else{
        auc <- AUC_of_cross_function(crossK) 
    }
    
    auc_sim <- c(auc_sim, auc)
    
    # summary of distance
    p_dist <- calculate_pairwise_distances_between_celltypes(
        image_real, c("islet", "endothelial"), "Cell.Type4")
    summary_p_dist <- calculate_summary_distances_between_celltypes(p_dist)
    if(!"islet/endothelial" %in% summary_p_dist$Pair) {
        avg_p_dist <- NA }else {
            avg_p_dist <- summary_p_dist[summary_p_dist$Pair  == "islet/endothelial", "Mean"]
        }
    
    m_dist <- calculate_minimum_distances_between_celltypes(
        image_real, cell_types_of_interest = c("islet", "endothelial"), feature_colname = "Cell.Type4")
    if(nrow(m_dist) == 0) {
        avg_m_dist1 <- NA
        avg_m_dist2 <- NA
    } else{
        summary_m_dist <- calculate_summary_distances_between_celltypes(m_dist)
        avg_m_dist1 <- summary_m_dist[summary_m_dist$Pair  == "islet/endothelial", "Mean"]
        avg_m_dist2 <- summary_m_dist[summary_m_dist$Pair  == "endothelial/islet", "Mean"]
    }
    
    avg_p_dist_real <- c(avg_p_dist_real, avg_p_dist)
    avg_m_dist1_real <- c(avg_m_dist1_real, avg_m_dist1)
    avg_m_dist2_real <- c(avg_m_dist2_real, avg_m_dist2)
    
    
    p_dist <- calculate_pairwise_distances_between_celltypes(
        image_sim, c("islet", "endothelial"), "Cell.Type2")
    summary_p_dist <- calculate_summary_distances_between_celltypes(p_dist)
    if(!"islet/endothelial" %in% summary_p_dist$Pair) {
        avg_p_dist <- NA }else {
            avg_p_dist <- summary_p_dist[summary_p_dist$Pair  == "islet/endothelial", "Mean"]
        }
    
    m_dist <- calculate_minimum_distances_between_celltypes(
        image_sim, cell_types_of_interest = c("islet", "endothelial"), feature_colname = "Cell.Type2")
    if(nrow(m_dist) == 0) {
        avg_m_dist1 <- NA
        avg_m_dist2 <- NA
    } else{
        summary_m_dist <- calculate_summary_distances_between_celltypes(m_dist)
        avg_m_dist1 <- summary_m_dist[summary_m_dist$Pair  == "islet/endothelial", "Mean"]
        avg_m_dist2 <- summary_m_dist[summary_m_dist$Pair  == "endothelial/islet", "Mean"]
    }
    
    avg_p_dist_sim <- c(avg_p_dist_sim, avg_p_dist)
    avg_m_dist1_sim <- c(avg_m_dist1_sim, avg_m_dist1)
    avg_m_dist2_sim <- c(avg_m_dist2_sim, avg_m_dist2)
    
    # CIN
    # the radius was computed within the loop when calculating ms and nms
    cin_real <- c(cin_real, average_percentage_of_cells_within_radius(image_real, reference_celltype = "islet",
                                              target_celltype = "endothelial", feature_colname = "Cell.Type4",
                                              radius = radius_real))

    cin_sim <- c(cin_sim, average_percentage_of_cells_within_radius(image_sim, reference_celltype = "islet",
                                                                      target_celltype = "endothelial", feature_colname = "Cell.Type2",
                                                                      radius = radius_sim))
    
}

final_metrics <- cbind(final_metrics, ms_real, ms_sim, nms_real, nms_sim, auc_real, 
                       auc_sim, avg_p_dist_real, avg_p_dist_sim, avg_m_dist1_real,
                       avg_m_dist2_real, avg_m_dist1_sim, avg_m_dist2_sim, cin_real, 
                       cin_sim)


# get the cell proportions of each image ####
p_beta_real <- c()
p_nonbeta_real <- c()
p_immune_real <- c()
p_endothelial_real <- c()
p_others_real <- c()

p_beta_sim <- c()
p_nonbeta_sim <- c()
p_immune_sim <- c()
p_endothelial_sim <- c()
p_others_sim <- c()

for (i in 1:845){
    image_real <- Real[[i]]
    image_sim <- even_sim_paired_final[[i]]
    
    p_real <- calculate_cell_proportions(image_real, plot.image = F,
                                         feature_colname = "Cell.Type")
    
    p_beta <- p_real[p_real$Cell_type == "beta", "Percentage"]
    if(length(p_beta) == 0) p_beta <- 0
    p_beta_real <- c(p_beta_real, p_beta)
    
    p_nonbeta <- sum(p_real[p_real$Cell_type %in% c("alpha", "gamma", "delta"), "Percentage"])
    if(length(p_nonbeta) == 0) p_nonbeta <- 0
    p_nonbeta_real <- c(p_nonbeta_real, p_nonbeta)
    
    p_immune <- sum(p_real[p_real$Cell_type %in% c("macrophage", "neutrophil", 
                                                   "otherimmune", "Tc", "Th", "B", "naiveTc"), 
                           "Percentage"])
    if(length(p_immune) == 0) p_immune <- 0
    p_immune_real <- c(p_immune_real, p_immune)
    
    p_endothelial <- p_real[p_real$Cell_type == "endothelial", "Percentage"]
    if(length(p_endothelial) == 0) p_endothelial <- 0
    p_endothelial_real <- c(p_endothelial_real, p_endothelial)
    
    p_others <- sum(p_real[p_real$Cell_type %in% c("acinar", "ductal", 
                                                   "stromal", "unknown"), 
                           "Percentage"])
    if(length(p_others) == 0) p_others <- 0
    p_others_real <- c(p_others_real, p_others)
    
    
    p_sim <- calculate_cell_proportions(image_sim, plot.image = F,
                                        feature_colname = "Cell.Type")
    p_beta <- p_sim[p_sim$Cell_type == "beta", "Percentage"]
    if(length(p_beta) == 0) p_beta <- 0
    p_beta_sim <- c(p_beta_sim, p_beta)
    
    p_nonbeta <- p_sim[p_sim$Cell_type == "nonbeta", "Percentage"]
    if(length(p_nonbeta) == 0) p_nonbeta <- 0
    p_nonbeta_sim <- c(p_nonbeta_sim, p_nonbeta)
    
    p_immune <- p_sim[p_sim$Cell_type == "immune", "Percentage"]
    if(length(p_immune) == 0) p_immune <- 0
    p_immune_sim <- c(p_immune_sim, p_immune)
    
    p_endothelial <- p_sim[p_sim$Cell_type == "endothelial", "Percentage"]
    if(length(p_endothelial) == 0) p_endothelial <- 0
    p_endothelial_sim <- c(p_endothelial_sim, p_endothelial)
    
    p_others <- p_sim[p_sim$Cell_type == "others", "Percentage"]
    if(length(p_others) == 0) p_others <- 0
    p_others_sim <- c(p_others_sim, p_others)
    
    print(i)
}

final_metrics <- cbind(final_metrics, p_beta_real, p_nonbeta_real, p_immune_real, 
                       p_endothelial_real, p_others_real, p_beta_sim, p_nonbeta_sim,
                       p_immune_sim, p_endothelial_sim, p_others_sim)
#####
save(final_metrics, file = "Objects/final_metrics.Rda")

