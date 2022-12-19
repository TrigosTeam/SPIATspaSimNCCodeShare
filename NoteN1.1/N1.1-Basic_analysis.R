library(SPIAT) #v1.0.4
# read file and format data ####
raw_inform_data <- "data/S6_[49209,17530]_cell_seg_data_updated.txt"

markers <- c("DAPI","CD3","PDL1","FOXP3","CD4","CD8","AMACR")

intensity_columns_interest <- c("Nucleus.DAPI..DAPI..Mean..Normalized.Counts..Total.Weighting.",
                                "Cytoplasm.CD3..Opal.520..Mean..Normalized.Counts..Total.Weighting.", 
                                "Membrane.PDL.1..Opal.540..Mean..Normalized.Counts..Total.Weighting.", 
                                "Cytoplasm.FOXP3..Opal.570..Mean..Normalized.Counts..Total.Weighting.", 
                                "Cytoplasm.CD4..Opal.620..Mean..Normalized.Counts..Total.Weighting.", 
                                "Cytoplasm.CD8..Opal.650..Mean..Normalized.Counts..Total.Weighting.", 
                                "Cytoplasm.AMACR..Opal.690..Mean..Normalized.Counts..Total.Weighting.")

#Formats an INFORM image into a SpatialExperiment class
#where the count assay stores the expression level of every marker (rows) for
#every cell (columns), and cell phenotype, x and y coordinates are stored
#under colData
formatted_image <- format_image_to_spe(format="inForm",
                                       path=raw_inform_data,
                                       markers=markers,
                                       dye_columns_interest=NULL,
                                       intensity_columns_interest=intensity_columns_interest)
# predict phenotypes ####
source("predict_phenotype_log_scale.R")
predicted_image <- predict_phenotype_log_scale(sce_object=formatted_image,
                                      thresholds = NULL,
                                      tumour_marker = "AMACR",
                                      baseline_markers = c("CD3", "CD4", "CD8"),
                                      reference_phenotypes=TRUE, nuclear_marker = "DAPI")


# use the phenotypes predicted by SPIAT
formatted_image <- predict_phenotypes(spe_object=formatted_image,
                                      thresholds = NULL,
                                      tumour_marker = "AMACR",
                                      baseline_markers = c("CD3", "CD4", "CD8"),
                                      reference_phenotypes=FALSE, nuclear_marker = "DAPI")

#Takes in the returned dataframe from marker_threshold_plot and
#generates a .pdf file containing scatter plots of actual expression and
#predicted expression for every marker.

# svglite("predicted_actual_plots.svg", width = 10, height = 4)
marker_prediction_plot(predicted_image, marker="CD8")

# dev.off()

# select and define cells ####

unique(formatted_image$Phenotype)

formatted_image <- select_celltypes(formatted_image, keep=TRUE,
                                    celltypes = c("AMACR", "CD3,CD4", "CD3,CD8",
                                                  "CD3,FOXP3,CD4", "PDL1", 
                                                  "AMACR,PDL1", "CD3"),
                                    feature_colname = "Phenotype")
formatted_image <- define_celltypes(formatted_image,
                                    categories = c("AMACR", "CD3,CD4", "CD3,CD8",
                                                   "CD3", "PDL1", "AMACR,PDL1", 
                                                   "CD3,FOXP3,CD4"),
                                    category_colname = "Phenotype",
                                    names = c("Tumour", "T_helper", "T_cyto", 
                                              "T_other", "PDL1", "Tumour_PDL1", 
                                              "Treg"),
                                    new_colname = "Cell.Type",print_names = T)

unique(formatted_image$Cell.Type)

##Convert to um from pixels
spatialCoords(formatted_image)[, "Cell.X.Position"] <- spatialCoords(formatted_image)[, "Cell.X.Position"]/2
spatialCoords(formatted_image)[, "Cell.Y.Position"] <- spatialCoords(formatted_image)[, "Cell.Y.Position"]/2

# cell proportions and counts ####

#Calculate the number and proportion of each cell phenotype in image
# svglite("cell_proportion.svg")
p_cells <- calculate_cell_proportions(spe_object = formatted_image, 
                                      feature_colname ="Cell.Type", 
                                      reference_celltypes = NULL)
# dev.off()

# plot the number of positive and negative cells for a marker
source("marker_boxplot_log_scale.R")
# svglite("boxplot.svg", height=2.5, width=3)
marker_boxplot_log_scale(formatted_image, "CD8")
# dev.off()

# plot cell categories ####
# svglite("plot_cell_category.svg", height=4.8, width=6)

# plot cells and markers ####
my_colors <- c("darkgrey", "skyblue", "darkcyan", "lightgreen", "orange")
plot_cell_categories(formatted_image, c("Tumour", "T_helper", "T_cyto", "T_other", "PDL1"), 
                     my_colors, "Cell.Type")
# dev.off()

# svglite("marker_level_heatmap.svg", height=4, width=5.5)
plot_marker_level_heatmap(formatted_image, num_splits = 100, "CD8")
# dev.off()

# svglite("cell_marker_level_heatmap.svg", height=4, width=5.5)
# plot_cell_marker_levels(formatted_image, marker="CD8")
# dev.off()

# Distances between phenotypes####
# svglite("distance_heatmap.svg", height=4, width=4.8)
minimum_distances <- calculate_minimum_distances_between_celltypes(
    formatted_image,cell_types_of_interest = c("Tumour","Treg", "T_helper", "T_cyto", "T_other", "PDL1"),
    feature_colname="Cell.Type")
summary_distance <- calculate_summary_distances_between_celltypes(minimum_distances)

#Takes the output of cell_distances and plot the distances as a heatmap
plot_distance_heatmap(summary_distance)
# dev.off()


##Calculate minimum distances (identifies the nearest cell type of the other type, and returns the distance)
distances_cyto_helper <- calculate_minimum_distances_between_celltypes(formatted_image,
                                                     cell_types_of_interest = c("T_cyto", "T_helper"),
                                                     feature_colname = "Cell.Type")

distances_cyto_PDL1 <- calculate_minimum_distances_between_celltypes(formatted_image,
                                                  cell_types_of_interest = c("T_cyto", "PDL1"),
                                                  feature_colname = "Cell.Type")

distances_cyto_helper <- distances_cyto_helper[distances_cyto_helper$RefType == "T_cyto" & 
                                                 distances_cyto_helper$NearestType == "T_helper",]

distances_cyto_PDL1 <- distances_cyto_PDL1[distances_cyto_PDL1$RefType == "T_cyto" & 
                                               distances_cyto_PDL1$NearestType == "PDL1",]

distances_cyto <- rbind(distances_cyto_helper, distances_cyto_PDL1)
distances_cyto$Pair <- paste(distances_cyto$RefType, distances_cyto$NearestType, sep="_")

# svglite("violin.svg", height=4, width=4.8)
ggplot(distances_cyto, aes(x=Pair, y=Distance))+
  geom_violin(aes(fill=Pair))+
  theme_bw()
# dev.off()

wilcox.test(distances_cyto$Dist[distances_cyto$Pair == "T_cyto_T_helper"],
            distances_cyto$Dist[distances_cyto$Pair == "T_cyto_PDL1"])$p.value

# tSNE ####
data <- data.frame(colData(formatted_image))
data$Cell.Type <- factor(data$Cell.Type , levels = c("T_cyto","PDL1","T_other", "T_helper","Tumour","Treg"))
colData(formatted_image)$Cell.Type <- data$Cell.Type

# svglite("tsne_no_other.svg")
try_image <- select_celltypes(formatted_image, c("PDL1","T_cyto","T_helper","Tumour","Treg"), feature_colname = "Cell.Type")

dimensionality_reduction_plot(try_image, plot_type = "TSNE", 
                              scale = T, feature_colname = "Cell.Type")
# dev.off()
