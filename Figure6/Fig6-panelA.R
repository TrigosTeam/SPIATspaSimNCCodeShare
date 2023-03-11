library(SPIAT)
load("../NoteN2/Objects/PCa_spe.Rda")

v1 <- c("AMACR" , "CD3,CD8" , "CD3,CD4" )
v2 <- c("Tumor", "Immune", "Immune")
v3 <- c("#FC8D62", "lightgray", "#7570B3") # colourblind friendly


spe <- PCa_spe[[13]][[1]]
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"


spe <- define_celltypes(spe, categories = v1, category_colname = "Phenotype", 
                          names = v2, new_colname = "Cell.Type")

svglite::svglite("Results/PanelA.svg", width = 5, height = 4)
g <- plot_cell_categories(spe, categories_of_interest = c( "Tumor", "Undefined", "Immune"), 
                          colour_vector = v3, 
                          feature_colname = "Cell.Type",
                          cex = 0.8, layered = T)
g
dev.off()
