load("../Figure4/Objects/Figure4_simulations.Rda")
library(SPIAT)
# gradient entropy #####
## Using Tumour cells as reference cells. (Panel B)
# input variables
gradient_pos <- seq(50, 800, 50)
cell_types_of_interest<-c("Tumour","Immune")
reference_marker = cell_types_of_interest[1]
target_marker = cell_types_of_interest[2]
feature_colname <- "Cell.Type"
categories_of_interest <- cell_types_of_interest
all_categories_of_interest <- c(categories_of_interest, "Others")

# output variables
entropy_gradient_list_all <- list()
peaks_all <- c()

# analysis 
for (i in 1:length(grad_imageL)){
    spe <- grad_imageL[[i]]
    # compute the entropy gradient and get the peak
    grad <- entropy_gradient_aggregated(spe_object = spe, 
                                        cell_types_of_interest = cell_types_of_interest,
                                        feature_colname = feature_colname, 
                                        radii = gradient_pos)
    entropy_gradient_list_all[[i]] <- grad
    
    v <- as.numeric(grad$gradient_df[1,3:(length(gradient_pos)+2)])
    
    # plot the gradient
    # svglite::svglite(paste("Results/entropy_", attr(spe, "name"),".svg", sep=""), width = 5, height = 4)
    plot(v, type = "b", lty = 2, pch = 16, cex = 1, ylim = c(0,1))
    # dev.off()
    # get the peak of the gradient
    peak <- grad$peak
    peaks_all[i] <- peak
}

# try mixing score #####

## Using Tumour cells as reference cells. (Panel D)

# input variables
gradient_pos <- seq(50, 800, 50)
cell_types_of_interest<-c("Tumour","Immune")
reference_marker = cell_types_of_interest[1]
target_marker = cell_types_of_interest[2]
feature_colname <- "Cell.Type"
categories_of_interest <- cell_types_of_interest
all_categories_of_interest <- c(categories_of_interest, "Others")

# output variables
nms_gradient_list_all <- list()

# define normalised mixing score function
nms <- function(spe_object, radius){
    ms_df <- mixing_score_summary(spe_object, 
                                  reference_celltype = reference_marker, 
                                  target_celltype = target_marker, 
                                  radius = radius,
                                  feature_colname = "Cell.Type")
    return(ms_df$Normalised_mixing_score)
}

# analysis
for (i in 1:length(grad_imageL)){
    spe <- grad_imageL[[i]]
    # compute the mixing score gradient
    grad <- compute_gradient(spe_object = spe, 
                             FUN = nms,
                             radii = gradient_pos)
    nms_gradient_list_all[[i]] <- unlist(grad)
    
    v <- unlist(grad)
    
    # plot the gradient
    # svglite::svglite(paste("Results/ms/", attr(spe, "name"),".svg", sep=""), width = 5, height = 4)
    plot(v, type = "b", lty = 2, pch = 16, cex = 1, ylim = c(0,2.2))
    # dev.off()
}

## Using Immune cells as reference cells. (Panel C)
# input variables
gradient_pos <- seq(50, 800, 50)
cell_types_of_interest<-c("Immune","Tumour")
reference_marker = cell_types_of_interest[1]
target_marker = cell_types_of_interest[2]
feature_colname <- "Cell.Type"
categories_of_interest <- cell_types_of_interest
all_categories_of_interest <- c(categories_of_interest, "Others")

# output variables
nms_gradient_list_all <- list()

# define normalised mixing score function
nms <- function(spe_object, radius){
    ms_df <- mixing_score_summary(spe_object, 
                                  reference_celltype = reference_marker, 
                                  target_celltype = target_marker, 
                                  radius = radius,
                                  feature_colname = "Cell.Type")
    return(ms_df$Normalised_mixing_score)
}

# analysis
for (i in 1:length(grad_imageL)){
    spe <- grad_imageL[[i]]
    # compute the mixing score gradient
    
    grad <- compute_gradient(spe_object = spe, 
                             FUN = nms,
                             radii = gradient_pos)
    nms_gradient_list_all[[i]] <- unlist(grad)
    
    v <- unlist(grad)
    
    # plot the gradient
    # svglite::svglite(paste("Results/nms_immune/", attr(spe, "name"),".svg", sep=""), width = 5, height = 4)
    plot(v, type = "b", lty = 2, pch = 16, cex = 1, ylim = c(0,3))
    # dev.off()
}
