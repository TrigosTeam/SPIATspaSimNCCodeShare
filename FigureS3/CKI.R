library(spaSim)
library(SPIAT)

# simulate scenarios when CKI gives a high value

CKI1 <- c()
for (i in 1:1000){
    image <- TIS(bg_sample = bg1, 
                 names_of_bg_cells = c("Others", "CellA", "CellB"),
                 proportions_of_bg_cells = c(0.7, 0.2, 0.1), 
                 plot_image = TRUE,
                 plot_colours = c("lightgray", "#D95F02", "#7570B3"))
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("CellA", "CellB"),
                                          feature_colname = "Cell.Type")
    
    CKI1 <- c(CKI1, crossing_of_crossK(df_cross))
    # save(CKI1, file = "Simulations/Objects/CKI1.Rda")
}

# 0.1006979 False positives
length(which(!is.na(CKI1)))/length(CKI1)
# sensitive to , shall we simulate a lot of this to calculate the probability of crossing 


# simulate few cells in the immune ring
CKI2 <- c()
for (i in 1:500){
    image <- TIS(
        bg_sample = bg1,
        n_immune_rings = 2,
        properties_of_immune_rings = list(
            I1 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Circle", 
                centre_loc = data.frame(x = 930, y = 1000), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.9)), 
            I2 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Oval",
                centre_loc = data.frame(x = 1330, y = 1100), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.9))),
        plot_image = TRUE,
        plot_categories = c("Cell_A", "Cell_B", "Others"),
        plot_colours = c( "#D95F02", "#7570B3", "lightgray")
    )
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("Cell_A", "Cell_B"),
                                          feature_colname = "Cell.Type")
    CKI2 <- c(CKI2, crossing_of_crossK(df_cross))
    # save(CKI2, file = "Objects/CKI2.Rda")
}

# Suppose the non-na values are false positives, then there are 0.2818004
length(which(!is.na(CKI2)))/length(CKI2)


# there are same proportions of B cells within the cluster and the immune ring
CKI2_2 <- c()
for (i in 1:500){
    image <- TIS(
        bg_sample = bg1,
        n_immune_rings = 2,
        properties_of_immune_rings = list(
            I1 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Circle", 
                centre_loc = data.frame(x = 930, y = 1000), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.9)), 
            I2 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Oval",
                centre_loc = data.frame(x = 1330, y = 1100), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.9))),
        plot_image = TRUE,
        plot_categories = c("Cell_A", "Cell_B", "Others"),
        plot_colours = c( "#D95F02", "#7570B3", "lightgray")
    )
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("Cell_A", "Cell_B"),
                                          feature_colname = "Cell.Type")
    CKI2_2 <- c(CKI2_2, crossing_of_crossK(df_cross))
    # save(CKI2_2, file = "Objects/CKI2_2.Rda")
}

# There is no false positives
length(which(!is.na(CKI2_2)))/length(CKI2_2)



# same Cell_B density everywhere
CKI2_3 <- c()
for (i in 1:500){
    image <- TIS(
        bg_sample = bg1,
        names_of_bg_cells = c("Cell_B", "Others"),
        proportions_of_bg_cells = c(0.1, 0.9),
        n_immune_rings = 2,
        properties_of_immune_rings = list(
            I1 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Circle", 
                centre_loc = data.frame(x = 930, y = 1000), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.9)), 
            I2 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Oval",
                centre_loc = data.frame(x = 1330, y = 1100), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.9))),
        plot_image = TRUE,
        plot_categories = c("Cell_A", "Cell_B", "Others"),
        plot_colours = c( "#D95F02", "#7570B3", "lightgray")
    )
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("Cell_A", "Cell_B"),
                                          feature_colname = "Cell.Type")
    CKI2_3 <- c(CKI2_3, crossing_of_crossK(df_cross))
    # save(CKI2_3, file = "Objects/CKI2_3.Rda")
}

# 0.188
length(which(!is.na(CKI2_3)))/length(CKI2_3)



# real Cell_B ring with no background
CKI3 <- c()
for (i in 1:500){
    image <- TIS(
        bg_sample = bg1,
        n_immune_rings = 2,
        properties_of_immune_rings = list(
            I1 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Circle", 
                centre_loc = data.frame(x = 930, y = 1000), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2)), 
            I2 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Oval",
                centre_loc = data.frame(x = 1330, y = 1100), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 20,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2))),
        plot_image = TRUE,
        plot_categories = c("Cell_A", "Cell_B", "Others"),
        plot_colours = c( "#D95F02", "#7570B3", "lightgray")
    )
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("Cell_A", "Cell_B"),
                                          feature_colname = "Cell.Type")
    CKI3 <- c(CKI3, crossing_of_crossK(df_cross))
    # save(CKI3, file = "Objects/CKI3.Rda")
}

# 0.786
length(which(!is.na(CKI3)))/length(CKI3)

# clear immune rings, multiple layers
CKI3_2 <- c()
for (i in 1:500){
    image <- TIS(
        bg_sample = bg1,
        n_immune_rings = 2,
        properties_of_immune_rings = list(
            I1 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Circle", 
                centre_loc = data.frame(x = 930, y = 1000), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 100,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2)), 
            I2 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Oval",
                centre_loc = data.frame(x = 1330, y = 1100), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 100,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2))),
        plot_image = TRUE,
        plot_categories = c("Cell_A", "Cell_B", "Others"),
        plot_colours = c( "#D95F02", "#7570B3", "lightgray")
    )
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("Cell_A", "Cell_B"),
                                          feature_colname = "Cell.Type")
    CKI3_2 <- c(CKI3_2, crossing_of_crossK(df_cross))
    # save(CKI3_2, file = "Objects/CKI3_2.Rda")
}

# 0.998
length(which(!is.na(CKI3_2)))/length(CKI3_2)


# multiple layer immune ring, with immune infiltration
CKI4 <- c()
for (i in 3:500){
    image <- TIS(
        bg_sample = bg1,
        n_immune_rings = 2,
        properties_of_immune_rings = list(
            I1 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Circle", 
                centre_loc = data.frame(x = 930, y = 1000), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 100,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2)), 
            I2 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Oval",
                centre_loc = data.frame(x = 1330, y = 1100), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 100,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2))),
        plot_image = TRUE,
        plot_categories = c("Cell_A", "Cell_B", "Others"),
        plot_colours = c( "#D95F02", "#7570B3", "lightgray")
    )
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("Cell_A", "Cell_B"),
                                          feature_colname = "Cell.Type")
    CKI4 <- c(CKI4, crossing_of_crossK(df_cross))
}
# save(CKI4, file = "Objects/CKI4.Rda")

# 0.914
length(which(!is.na(CKI4)))/length(CKI4)



# multiple layer immune ring, with immune cells everywhere
CKI4_2 <- c()
for (i in 1:500){
    image <- TIS(
        bg_sample = bg1,
        names_of_bg_cells = c("Cell_B", "Others"),
        proportions_of_bg_cells = c(0.1, 0.9),
        n_immune_rings = 2,
        properties_of_immune_rings = list(
            I1 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Circle", 
                centre_loc = data.frame(x = 930, y = 1000), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 100,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2)), 
            I2 = list(
                name_of_cluster_cell = "Cell_A", size = 400, shape = "Oval",
                centre_loc = data.frame(x = 1330, y = 1100), 
                infiltration_types = c("Cell_B", "Others"), 
                infiltration_proportions = c(0.1, 0.05),
                name_of_ring_cell = "Cell_B", immune_ring_width = 100,
                immune_ring_infiltration_types = c("Others"), 
                immune_ring_infiltration_proportions = c(0.2))),
        plot_image = TRUE,
        plot_categories = c("Cell_A", "Cell_B", "Others"),
        plot_colours = c( "#D95F02", "#7570B3", "lightgray")
    )
    
    df_cross <- calculate_cross_functions(image, 
                                          cell_types_of_interest = c("Cell_A", "Cell_B"),
                                          feature_colname = "Cell.Type")
    CKI4_2 <- c(CKI4_2, crossing_of_crossK(df_cross))
}
# save(CKI4_2, file = "Objects/CKI4_2.Rda")

# 0.994
length(which(!is.na(CKI4_2)))/length(CKI4_2)
