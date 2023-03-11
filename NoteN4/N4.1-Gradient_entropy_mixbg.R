library(spaSim)
library(SPIAT)
library(ggplot2)
library(svglite)

# Simulate images with fixed Cell_type1 cells and increasing Cell_type2 cells ####
background_sample <- bg1
names_of_cell_types <- c("Cell_type1" , "Cell_type2", "Others")
props_type1 <- rep(0.4, 11) # fixed proportion of Cell_type1 cells 0.4
props_type2 <- seq(0.05, 0.55, 0.05) # increasing proportions of Cell_type2 cells
props_Others <- seq(0.55, 0.05, -0.05) # decreasing proportions of Other cells
proportions_of_cell_types <- list(props_type1,props_type2, props_Others)
# pdf("~/Documents/paper/SPIAT_figs/S1/sim_images.pdf")
set.seed(610)
imageL <- multiple_background_images(background_sample,
                                     names_of_cell_types,
                                     proportions_of_cell_types,
                                     plot_image = TRUE,
                                     plot_colours = c("#D95F02", "#7570B3", "lightgray"))
# dev.off()

# entropy_gradient ####
# use Cell_type1 as reference
# input variables
gradient_pos <- seq(50, 800, 50)
cell_types_of_interest<-c("Cell_type1","Cell_type2")
reference_marker = cell_types_of_interest[1]
target_marker = cell_types_of_interest[2]
feature_colname <- "Cell.Type"
categories_of_interest <- cell_types_of_interest
all_categories_of_interest <- c(categories_of_interest, "Others")

# output variables
entropy_gradient_list_all <- list()
peaks_all <- c()


for (i in 1:length(imageL)){
    sce <- imageL[[i]]
    attr(sce, "name") <- i
    # compute the entropy gradient and get the peak
    grad <- entropy_gradient_aggregated(spe_object = sce, 
                                        cell_types_of_interest = cell_types_of_interest,
                                        feature_colname = feature_colname, 
                                        radii = gradient_pos)
    entropy_gradient_list_all[[i]] <- grad
    
    v <- as.numeric(grad$gradient_df[1,3:(length(gradient_pos)+2)])
    
    # plot the gradient
    # svglite::svglite(paste("Results/more_simulations_Type1Ref/entropy_", attr(sce, "name"),".svg", sep=""), width = 5, height = 4)
    plot(v, type = "b", lty = 2, pch = 16, cex = 1, ylim = c(0,1))
    # dev.off()
    # get the peak of the gradient
    peak <- grad$peak
    peaks_all[i] <- peak
}

# use Cell_type1 as reference
# input variables
gradient_pos <- seq(50, 800, 50)
cell_types_of_interest<-c("Cell_type2","Cell_type1")
reference_marker = cell_types_of_interest[1]
target_marker = cell_types_of_interest[2]
feature_colname <- "Cell.Type"
categories_of_interest <- cell_types_of_interest
all_categories_of_interest <- c(categories_of_interest, "Others")

# output variables
entropy_gradient_list_all <- list()
peaks_all <- c()


for (i in 1:length(imageL)){
    sce <- imageL[[i]]
    attr(sce, "name") <- i
    # compute the entropy gradient and get the peak
    grad <- entropy_gradient_aggregated(spe_object = sce, 
                                        cell_types_of_interest = cell_types_of_interest,
                                        feature_colname = feature_colname, 
                                        radii = gradient_pos)
    entropy_gradient_list_all[[i]] <- grad
    
    v <- as.numeric(grad$gradient_df[1,3:(length(gradient_pos)+2)])
    
    # plot the gradient
    # svglite::svglite(paste("Results/more_simulations_Type2Ref/entropy_", attr(sce, "name"),".svg", sep=""), width = 5, height = 4)
    plot(v, type = "b", lty = 2, pch = 16, cex = 1, ylim = c(0,1))
    # dev.off()
    # get the peak of the gradient
    peak <- grad$peak
    peaks_all[i] <- peak
}


for (i in 1:length(imageL)){
    sce <- imageL[[i]]
    # compute the entropy gradient and get the peak
    print(table(sce$Cell.Type))
}

