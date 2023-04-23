# This script is for the scatter plots of the metrics
library(ggplot2)
load("Objects/final_metrics.Rda")

metrics <- c("width", "height", "n_cluster_cells", "R_BC", "ms", "nms", "auc",
             "avg_p_dist", "avg_m_dist1", "avg_m_dist2", "cin", "p_beta", "p_nonbeta",
             "p_immune", "p_endothelial", "p_others")


for (metric in metrics){
    metric_real <- paste0(metric, "_real")
    metric_sim <- paste0(metric, "_sim")
    
    fit <- lm( final_metrics[, metric_real] ~ final_metrics[, metric_sim])
    print(metric)
    print(summary(fit))
    # svg(paste0("Results/", metric, ".svg"), width = 3, height = 3)
    g <- ggplot(final_metrics, aes(get(metric_real), get(metric_sim))) +
        geom_point() +
        geom_smooth(method = "lm") +
        ggtitle(metric)
    plot(g)
    # dev.off()
}

metrics_log <- c( "ms", "nms")

for (metric in metrics_log){
    metric_real <- paste0(metric, "_real")
    metric_sim <- paste0(metric, "_sim")
    
    fit <- lm( final_metrics[, metric_real] ~ final_metrics[, metric_sim])
    print(metric)
    print(summary(fit))
    # svg(paste0("Results/", metric, "_log.svg"), width = 3, height = 3)
    g <- ggplot(final_metrics, aes(log(get(metric_real)), log(get(metric_sim)))) +
        geom_point() +
        geom_smooth(method = "lm") +
        ggtitle(metric)
    plot(g)
    # dev.off()
}





# R-square
# avg_p_dist: 0.90
# avg_min_dist_1: 0.33
# avg_min_dist_2: 0.08
# CIN: 0.52
# auc: 0.24
# ms: 0.47
# nms: 0.49

## Plot the paired real and simulated images

load("Objects/diabetes_real_spe.Rda")
load("Objects/even_sim_paired_final.Rda")
spe <- even_sim_paired_final[[3]]
library(SPIAT)
plot_cell_categories(spe, categories_of_interest = c("beta", "nonbeta", "immune", "endothelial", "others"),
                     colour_vector = c("#E7298A", "#D95F02", "#1B9E77", "#7570B3", "lightgray" ),
                     feature_colname = "Cell.Type")

spe2 <- Real[[3]]

spe2 <- define_celltypes(spe2, categories = c("beta", "nonbeta", "ductal", "stromal", 
                                              "unknown", "acinar", "endothelial", "Tc", "macrophage", 
                                              "neutrophil", "otherimmune", "B"), category_colname = "Cell.Type3",
                         names = c("beta", "nonbeta", "others", "others", "others", 
                                   "others", "endothelial", "immune", "immune", "immune", "immune", "immune"),
                         new_colname = "Cell.Type5")
plot_cell_categories(spe2, categories_of_interest = c("beta", "nonbeta", "immune", "endothelial", "others"),
                     colour_vector = c("#E7298A", "#D95F02", "#1B9E77", "#7570B3", "lightgray" ),
                     feature_colname = "Cell.Type5")
