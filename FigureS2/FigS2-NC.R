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
# pdf("Results/sim_images.pdf")
set.seed(610)
imageL <- multiple_background_images(background_sample,
                                     names_of_cell_types,
                                     proportions_of_cell_types,
                                     plot_image = TRUE,
                                     plot_colours = c("#D95F02", "#7570B3", "lightgray"))
# dev.off()
# Compute mixing score and normalised mixing score ####
## use Immune cells as target
reference_celltype <- "Cell_type1"
target_celltype <- "Cell_type2"
feature_colname <- "Cell.Type"
ms_df <- data.frame()
for (image in imageL) {
    out <- mixing_score_summary(image, reference_celltype, target_celltype, 
                                  feature_colname, radius = 20)
    ms_df <- rbind(ms_df, out)
}
row.names(ms_df) <- 1:11

df <- reshape2::melt(ms_df, id.vars = c("Reference", "Target","Number_of_target_cells",
                                        "Number_of_reference_cells", "Reference_target_interaction",
                                        "Reference_reference_interaction") , 
                     value.vars = c("Mixing_score","Normalised_mixing_score"),
                     value.name = "value")

svglite("Results/(N)MS-target.svg", width = 5, height = 4)
ggplot(df, aes(x = Number_of_target_cells, y = value)) +
    geom_line(linetype = "longdash", colour = "red")+
    geom_point() +
    facet_grid(variable ~ .,scales="fixed") +
    theme(axis.text.x = element_text(angle = 0)) 
dev.off()

## use Immune cells as reference
reference_celltype <- "Cell_type2"
target_celltype <- "Cell_type1"
ms_df2 <- data.frame()
for (image in imageL) {
    out <- mixing_score_summary(image, reference_celltype, target_celltype, 
                                feature_colname, radius = 20)
    ms_df2 <- rbind(ms_df2, out)
}
row.names(ms_df2) <- 1:11
df2 <- reshape2::melt(ms_df2, id.vars = c("Reference", "Target","Number_of_target_cells",
                                        "Number_of_reference_cells", "Reference_target_interaction",
                                        "Reference_reference_interaction") , 
                     value.vars = c("Mixing_score","Normalised_mixing_score"),
                     value.name = "value")
svglite("Results/(N)MS-reference.svg", width = 5, height = 4)
ggplot(df2, aes(x = Number_of_reference_cells, y = value)) +
    geom_line(linetype = "longdash", colour = "red")+
    geom_point() +
    facet_grid(variable ~ .,scales="fixed") +
    theme(axis.text.x = element_text(angle = 0)) 
dev.off()



