library(spaSim)
library(SPIAT)
library(ggplot2)
# simulate infiltration ####
set.seed(1995)
image <- multiple_background_images(bg_sample = bg1,
                                           props = list(0.1,0.1,0.8))
image <- image[[1]]

set.seed(995)
images <- multiple_images_with_immune_rings(bg_sample = image,
                                                  cluster_size = 200,
                                                  ring_shape = 2,
                                                  prop_infiltration = c(0.2, 0.3, 0.4, 0.5),
                                                  ring_width = 80,
                                                  cluster_loc_x = 0,
                                                  cluster_loc_y = 0,
                                                  prop_ring_infiltration = 0, 
                                            plot_image = T,
                                            plot_categories = c("Tumour", "Immune", "Others"),
                                            plot_colours = c("#D95F02", "#7570B3", "lightgray"))

df_list <- list()
a <- 0
for (spe in images){
    spe <-identify_bordering_cells(spe, "Tumour", "Cell.Type", n_to_exclude = 30)
    spe <- calculate_distance_to_margin(spe)
    spe <- define_structure(spe, "Immune")
    a <- a+1
    df_list[[a]] <- calculate_proportions_of_cells_in_structure(spe, "Immune", "Cell.Type")
}

metric <- c()
for (i in 1:4){
    metric <- c(metric, df_list[[i]][["P.Infiltrated.CoI"]][1])
}
data <- data.frame(x = 1:4, metric)

ggplot(data, aes(x, metric)) + 
    geom_point()+
    geom_line(linetype = "longdash", colour = "red")


# simulate separation ####
set.seed(1995)
image <- multiple_background_images(bg_sample = bg1,
                                    props = list(0.1,0.1,0.8))
image <- image[[1]]
set.seed(0)
images <- multiple_images_with_clusters(bg_sample = image, 
                                            cluster_shape = 2, 
                                              prop_infiltration = 0.1, 
                                              cluster_size = 200, 
                                              cluster_loc_x = -400, 
                                              cluster_loc_y = -350)
image <- images[[1]]
set.seed(0)
images <- multiple_images_with_clusters(bg_sample = image,
                                              cluster_shape = 3, 
                                              prop_infiltration = 0, 
                                              cluster_size = 400, 
                                              cluster_loc_x = c(-40, 50, 200, 320),
                                              cluster_loc_y = c(-40, 100, 200, 300),
                                        plot_image = T,
                                        plot_categories = c("Tumour", "Immune", "Others"),
                                        plot_colours = c("#D95F02", "#7570B3", "lightgray"))

df_list <- list()
a <- 0
for (spe in images){
    spe <-identify_bordering_cells(spe, "Tumour", "Cell.Type", n_to_exclude = 30, ahull_alpha = 50)
    spe <- calculate_distance_to_margin(spe)
    spe <- define_structure(spe, "Immune")
    a <- a+1
    df_list[[a]] <- calculate_summary_distances_of_cells_to_borders(spe, "Immune", "Cell.Type")
}


metric <- c()
for (i in 1:4){
    metric <- c(metric, df_list[[i]][["Mean_d"]][4])
}
data <- data.frame(x = 1:4, metric)
ggplot(data, aes(x, metric)) + 
    geom_point()+
    geom_line(linetype = "longdash", colour = "red")
