library(SPIAT)
load("../NoteN2/Objects/PCa_spe.Rda")

v1 <- c("AMACR" , "CD3,CD8" , "CD3,CD4" )
v2 <- c("Tumor", "Immune", "Immune")
v3 <- c("red", "darkgreen", "lightgray")
v4 <- c("Tumor", "Immune", "Undefined")

spe_1 <- PCa_spe[[6]][[6]]
spe_1$Phenotype[which(is.na(spe_1$Phenotype))] <- "OTHER"
spe_2 <- PCa_spe[[13]][[8]]
spe_2$Phenotype[which(is.na(spe_2$Phenotype))] <- "OTHER"
spe_3 <- PCa_spe[[2]][[4]]
spe_3$Phenotype[which(is.na(spe_3$Phenotype))] <- "OTHER"

spe_1 <- define_celltypes(spe_1, categories = v1, category_colname = "Phenotype", 
                          names = v2, new_colname = "Cell.Type")
spe_2 <- define_celltypes(spe_2, categories = v1, "Phenotype", names = v2, "Cell.Type")
spe_3 <- define_celltypes(spe_3, categories = v1, "Phenotype", names = v2, "Cell.Type")

svglite::svglite("Results/real_imgs/clusters.svg", width = 5, height = 4)
plot_cell_categories(spe_1, categories_of_interest = c( "Tumor", "Undefined", "Immune"), 
                     colour_vector = c("red", "lightgray", "darkgreen"), 
                     feature_colname = "Cell.Type",
                     cex = 0.8, layered = T)
dev.off()


svglite::svglite("Results/real_imgs/rings.svg", width = 5, height = 4)
plot_cell_categories(spe_2, categories_of_interest = c( "Tumor", "Undefined", "Immune"), 
                     colour_vector = c("red", "lightgray", "darkgreen"), 
                     feature_colname = "Cell.Type",
                     cex = 0.8, layered = T)
dev.off()


svglite::svglite("Results/real_imgs/vessels.svg", width = 5, height = 4)
plot_cell_categories(spe_3, categories_of_interest = c( "Tumor", "Undefined", "Immune"), 
                     colour_vector = c("red", "lightgray", "darkgreen"), 
                     feature_colname = "Cell.Type",
                     cex = 0.8, layered = T)
dev.off()
