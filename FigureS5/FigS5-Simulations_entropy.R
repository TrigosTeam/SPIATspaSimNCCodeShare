library(spaSim)
library(SPIAT)

# Simulations
prop_red <- c(0.7, 0.5, 0.3, 0.2, 0.1)
ImageList <- list()
set.seed(610)
for (i in 1:5){
    p <- prop_red[i]
    image <- TIS(bg_sample = bg1, 
                 names_of_bg_cells = c("Cell A", "Cell B", "Others"),
                 proportions_of_bg_cells = c(p, 0.3, 0.7-p), 
                 plot_image = TRUE, 
                 plot_categories = c("Cell A", "Cell B", "Others"),
                 plot_colours = c("#D95F02", "#7570B3", "lightgray"))
    ImageList[[i]] <- image
    # save(ImageList, file = "Objects/S5_Simulations.Rda")
}

# Grid entropy
gridList <- list()
for (i in 1:5){
    image <- ImageList[[i]]
    grid <- grid_metrics(spe_object = image, FUN = calculate_entropy, n_split = 20,
                         cell_types_of_interest = c("Cell A", "Cell B"), 
                 feature_colname = "Cell.Type")
    gridList[[i]] <- grid
}

# Plot density
for (i in 1:5){
    grid <- gridList[[i]]
    plot(density(grid@data@values, na.rm = T), ylim = c(0,8))
}
