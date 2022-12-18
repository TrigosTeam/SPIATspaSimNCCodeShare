## This script is for reading in diabetes data and extract metrics from each image.
library(SPIAT)
library(ggplot2)
library("RColorBrewer")
library(dplyr)
# Read image ####

# Publically available data downloaded from https://data.mendeley.com/datasets/cydmwsfztj/2
# Data was read with `read.delim()` from public dataset and saved as `.Rdata` for further analysis.

#data <- read.delim("/home/atrigos/Spatial_paper/Public_datasets/Diabetes_mass_cyto/All_Cells.csv",sep=",")
#save(data, file="/home/atrigos/Spatial_paper/Public_datasets/Diabetes_mass_cyto/data.Rdata")
load("data/data.Rdata")

metadata <- read.delim("data/Version_1/Metadata.csv",sep=",")
panel <- read.delim("data/Version_1/Panel.csv",sep=",")
celltypes <- read.delim("data/Version_2/CellTypes.csv",sep=",")
donors <- read.delim("data/Version_2/Donors.csv",sep=",")

celltypes$cell <- substr(celltypes$id, 5, nchar(celltypes$id))
celltypes$region <- substr(celltypes$id, 1, 1)
celltypes$patient <- donors$case[match(celltypes$region, donors$slide)]

image_map <- unique(celltypes[,c("core", "region", "patient")])
image_map$image <- 1:nrow(image_map)
image_map$Image_ID <- paste(image_map$patient,
                            image_map$region,
                            image_map$image, sep="_")

celltypes$Image_ID <- image_map$Image_ID[match(celltypes$core, image_map$core)]
celltypes$CellType2 <- celltypes$CellType
celltypes$CellType2[celltypes$CellType2 %in% c("alpha", "delta", "beta",  "gamma")] <- "islet"
data$Image_ID <- image_map$Image_ID[match(data$ImageNumber, image_map$image)]

# Extract metrics####
#i <- "6362_A_1"
for(i in unique(data$Image_ID)){
  n_cells <- list()
  image_size <- list()
  proportion_cells <- list()
  min_distance <- list()
  overall_min_distance <- list()
  cin <- list()
  mixing <- list()
  entropy <- list()
  cell_proportions_in_structure <- list()
  cell_proportions_in_structure_islet_cells <- list()
  cell_proportions_in_islet <- list()
  distances_to_border <- list()
  
  temp_image <- data[data$Image_ID == i,]
  temp_celltypes <- celltypes[celltypes$Image_ID == i,]
  
  temp_image$CellType <- temp_celltypes$CellType[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  temp_image$CellType2 <- temp_celltypes$CellType2[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  temp_image$CellCat <- temp_celltypes$CellCat[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  
  map <- data.frame(CellType = temp_image$CellType,
                    CellType2 = temp_image$CellType2,
                    CellCat = temp_image$CellCat)
  map <- unique(map)
  
  formatted_image <- format_image_to_spe(format = "general", 
                                         phenotypes = temp_image$CellType, 
                                         coord_x = temp_image$Location_Center_X,
                                         coord_y = temp_image$Location_Center_Y)
  rownames(colData(formatted_image)) <- paste("Cell", 1:nrow(colData(formatted_image)), sep="_")
  
  
  formatted_image <- define_celltypes(
    formatted_image, 
    categories = map$CellType, 
    category_colname = "Phenotype", 
    names = map$CellCat,
    new_colname = "Cell.Cat")
  
  formatted_image <- define_celltypes(
    formatted_image, 
    categories = map$CellType, 
    category_colname = "Phenotype", 
    names = map$CellType,
    new_colname = "Cell.Type")
  
  formatted_image <- define_celltypes(
    formatted_image, 
    categories = map$CellType, 
    category_colname = "Phenotype", 
    names = map$CellType2,
    new_colname = "Cell.Type2")
  
  #Total number of cells
  n_cells[[i]] <- ncol(formatted_image)
  
  #Image size
  image_size[[i]] <- c(X_coord = max(spatialCoords(formatted_image)[,1]), Y_coord = max(spatialCoords(formatted_image)[,2]))
  
  ##Proportion of cells in the entire image. You can get islet size from here (see my code for the paper below)
  proportion_cells[[i]] <- calculate_cell_proportions(formatted_image, 
                                                      reference_celltypes=NULL, 
                                                      feature_colname ="Cell.Type",
                                                      plot.image = FALSE)
  
  #Minimum distance between cell pairs - for later
  min_dist <- calculate_minimum_distances_between_celltypes(
    formatted_image, 
    cell_types_of_interest = unique(colData(formatted_image)$Cell.Type),
    feature_colname = "Cell.Type")
  
  min_distance[[i]] <- calculate_summary_distances_between_celltypes(min_dist)
  
  #Radius for various calculations
  radius <- average_minimum_distance(formatted_image)*3
  
  overall_min_distance[[i]] <- average_minimum_distance(formatted_image)
  
  cells <- unique(colData(formatted_image)$Cell.Type)
  cells <- cells[!(cells %in% c("unknown", "otherimmune"))]
  cells <- c(cells, "islet")

  ##Identify margin cells
  formatted_border <- identify_bordering_cells(formatted_image, 
                                                reference_cell = "islet", ahull_alpha = 20,
                                                feature_colname="Cell.Type2", 
                                                n_of_polygons = 1, draw=FALSE)

  formatted_distance <- calculate_distance_to_margin(formatted_border)
  
  ##Cells in structure
  non_islet <- c("acinar","ductal","stromal","endothelial","Th","macrophage","neutrophil",
                 "Tc","naiveTc","B")
  islet <- c("alpha", "delta", "beta", "gamma")
  
  ##Non-islet cells in structures
  formatted_structure <- define_structure(
    formatted_distance, cell_types_of_interest = non_islet, 
    feature_colname = "Cell.Type", n_margin_layers = 5)
  
  cell_proportions_in_structure[[i]] <- calculate_proportions_of_cells_in_structure(
    formatted_structure, cell_types_of_interest = non_islet, 
    feature_colname ="Cell.Type")

  ##Composition of islets
  formatted_structure <- define_structure(
    formatted_distance, cell_types_of_interest = islet, 
    feature_colname = "Cell.Type", n_margin_layers = 5)
  
  cell_proportions_in_islet[[i]] <- table(formatted_structure@colData[,c("Cell.Type", "Region")])
  
  #Distance to the margin (all cell types)
  names_of_other_cells <- unique(colData(formatted_image)$Cell.Type)
  
  distances_to_border[[i]] <- calculate_summary_distances_of_cells_to_borders(
    formatted_structure, cell_types_of_interest = names_of_other_cells, 
    "Cell.Type")
  
  # save(n_cells, file=paste("~/Documents/paper/paper_new_figs/Objects/n_cells/n_cells_", i, ".Rdata", sep=""))  
  # save(image_size, file=paste("~/Documents/paper/paper_new_figs/Objects/image_size/image_size_", i, ".Rdata", sep=""))  
  # save(proportion_cells, file=paste("~/Documents/paper/paper_new_figs/Objects/proportion_cells/proportion_cells_", i, ".Rdata", sep=""))
  # save(min_distance, file=paste("~/Documents/paper/paper_new_figs/Objects/min_distance/min_distance_", i, ".Rdata", sep=""))
  # save(overall_min_distance, file=paste("~/Documents/paper/paper_new_figs/Objects/overall_min_distance/overall_min_distance_", i, ".Rdata", sep=""))
  # save(cell_proportions_in_structure, file=paste("~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_structure/cell_proportions_in_structure_", i, ".Rdata", sep=""))
  # save(cell_proportions_in_islet, file=paste("~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_islet/cell_proportions_in_islet_", i, ".Rdata", sep=""))
  # save(distances_to_border, file=paste("~/Documents/paper/paper_new_figs/Objects/distances_to_border/distances_to_border_", i, ".Rdata", sep=""))
  print(i)
}


# Collapsing metrics to single data frames ####
n_cells_all <- vector()
image_size_all <- vector()
proportion_cells_all <- vector()
overall_min_distance_all <- vector()
min_distance_all <- vector()
cell_proportions_in_structure_all <- vector()
distances_to_border_all <- vector()
cell_proportions_in_islet_all <- vector()
for(image_ID in unique(data$Image_ID)){
  
  patient <- strsplit(image_ID, "_")[[1]][1]
  slide <- strsplit(image_ID, "_")[[1]][2]
  image <- strsplit(image_ID, "_")[[1]][3]
  
  
  # load(paste("~/Documents/paper/paper_new_figs/Objects/n_cells/n_cells_", image_ID, ".Rdata", sep=""))  
  n_cells_all <- as.data.frame(rbind(n_cells_all,
                       c(Total_cells = n_cells[[1]], Patient = patient, Slide = slide, Image = image)))
  n_cells_all$Total_cells <- as.numeric(as.character(n_cells_all$Total_cells))
  
  # load(paste("~/Documents/paper/paper_new_figs/Objects/image_size/image_size_", image_ID, ".Rdata", sep=""))  
  image_size_all <- as.data.frame(rbind(image_size_all,
                          c(image_size[[1]], Patient = patient, Slide = slide, Image = image)))
  image_size_all$X_coord <- as.numeric(as.character(image_size_all$X_coord))
  image_size_all$Y_coord <- as.numeric(as.character(image_size_all$Y_coord))
  
  # load(paste("~/Documents/paper/paper_new_figs/Objects/overall_min_distance/overall_min_distance_", image_ID, ".Rdata", sep=""))
  overall_min_distance_all <- as.data.frame(rbind(overall_min_distance_all,
                                    c(Overall_min_distance=overall_min_distance[[1]], Patient = patient, Slide = slide, Image = image)))
  overall_min_distance_all$Overall_min_distance <- as.numeric(as.character(overall_min_distance_all$Overall_min_distance))
  
  # load(paste("~/Documents/paper/paper_new_figs/Objects/proportion_cells/proportion_cells_", image_ID, ".Rdata", sep=""))
    proportion_cells_all <- rbind(proportion_cells_all,
                                cbind(proportion_cells[[1]], Patient = patient, Slide = slide, Image = image))
  
  # load(paste("~/Documents/paper/paper_new_figs/Objects/min_distance/min_distance_", image_ID, ".Rdata", sep=""))
  min_distance_all <- rbind(min_distance_all,
                            cbind(min_distance[[1]], Patient = patient, Slide = slide, Image = image))
  

  # load(paste("~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_structure/cell_proportions_in_structure_", image_ID, ".Rdata", sep=""))
  cell_proportions_in_structure_all <- rbind(cell_proportions_in_structure_all,
                                             cbind(cell_proportions_in_structure[[1]], Patient = patient, Slide = slide, Image = image))
  
  # load(paste("~/Documents/paper/paper_new_figs/Objects/distances_to_border/distances_to_border_", image_ID, ".Rdata", sep=""))
  distances_to_border_all <- rbind(distances_to_border_all,
                                   cbind(distances_to_border[[1]], Patient = patient, Slide = slide, Image = image))
  
  # load(paste("~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_islet/cell_proportions_in_islet_", image_ID, ".Rdata", sep=""))
  
  if(length(colnames(cell_proportions_in_islet[[1]])) == 3){
    cell_proportions_in_islet_all <- as.data.frame(rbind(cell_proportions_in_islet_all,
                                           cbind(Cell.Type = rownames(cell_proportions_in_islet[[1]]), cell_proportions_in_islet[[1]], Patient = patient, Slide = slide, Image = image)))
  }else{
    temp <- cell_proportions_in_islet[[1]]
    NAs <- rep(NA, nrow(temp))
    if(!("Outside" %in% colnames(temp))){
      cell_proportions_in_islet[[1]] <- cbind(cell_proportions_in_islet[[1]], Outside=NAs)
    }
    if(!("Inside" %in% colnames(temp))){
      cell_proportions_in_islet[[1]] <- cbind(cell_proportions_in_islet[[1]], Inside=NAs)
    }
    if(!("Border" %in% colnames(temp))){
      cell_proportions_in_islet[[1]] <- cbind(cell_proportions_in_islet[[1]], Border=NAs)
    }
    cell_proportions_in_islet_all <- as.data.frame(rbind(cell_proportions_in_islet_all,
                                           cbind(Cell.Type = rownames(cell_proportions_in_islet[[1]]), cell_proportions_in_islet[[1]], Patient = patient, Slide = slide, Image = image)))
    
  }
  cell_proportions_in_islet_all$Border <- as.numeric(as.character(cell_proportions_in_islet_all$Border))
  cell_proportions_in_islet_all$Inside <- as.numeric(as.character(cell_proportions_in_islet_all$Inside))
  cell_proportions_in_islet_all$Outside <- as.numeric(as.character(cell_proportions_in_islet_all$Outside))
  
  print(image)
}

save(n_cells_all, file="~/Documents/paper/paper_new_figs/Objects/n_cells_all.Rdata")  
save(image_size_all, file="~/Documents/paper/paper_new_figs/Objects/image_size_all.Rdata")  
save(proportion_cells_all, file="~/Documents/paper/paper_new_figs/Objects/proportion_cells_all.Rdata")
save(overall_min_distance_all, file="~/Documents/paper/paper_new_figs/Objects/overall_min_distance_all.Rdata")
save(min_distance_all, file="~/Documents/paper/paper_new_figs/Objects/min_distance_all.Rdata")
save(cell_proportions_in_structure_all, file="~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_structure_all.Rdata")
save(distances_to_border_all, file="~/Documents/paper/paper_new_figs/Objects/distances_to_border_all.Rdata")
save(cell_proportions_in_islet_all, file="~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_islet_all.Rdata")


# ANALYSES FOR PAPER ####
## Non-islet cells in structures ####
load("~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_structure_all.Rdata")

cell_proportions_in_structure_all$Stage <- donors$stage[match(cell_proportions_in_structure_all$Slide, donors$slide)]

cell_proportions_all_cells_structure <- cell_proportions_in_structure_all[cell_proportions_in_structure_all$Relative_to == "All_cells_in_the_structure",]
cell_proportions_across_images <- cell_proportions_in_structure_all[cell_proportions_in_structure_all$Relative_to == "The_same_cell_type_in_the_whole_image",]
  
cell_proportions_all_cells_structure_agg_p_infiltrated <- aggregate(P.Infiltrated.CoI ~ Patient+Stage+Cell.Type, cell_proportions_all_cells_structure, mean)
cell_proportions_all_cells_structure_agg_p_internal <- aggregate(P.Internal.Margin.CoI ~ Patient+Stage+Cell.Type, cell_proportions_all_cells_structure, mean)
cell_proportions_all_cells_structure_agg_p_external <- aggregate(P.External.Margin.CoI ~ Patient+Stage+Cell.Type, cell_proportions_all_cells_structure, mean)
cell_proportions_all_cells_structure_agg_p_stromal <- aggregate(P.Stromal.CoI ~ Patient+Stage+Cell.Type, cell_proportions_all_cells_structure, mean)

cell_proportions_all_cells_structure_agg_p_infiltrated$Location <- "Infiltrated"
cell_proportions_all_cells_structure_agg_p_internal$Location <- "Internal"
cell_proportions_all_cells_structure_agg_p_external$Location <- "External"
cell_proportions_all_cells_structure_agg_p_stromal$Location <- "Stromal"

colnames(cell_proportions_all_cells_structure_agg_p_infiltrated)[4] <- "Proportion"
colnames(cell_proportions_all_cells_structure_agg_p_internal)[4] <- "Proportion"
colnames(cell_proportions_all_cells_structure_agg_p_external)[4] <- "Proportion"
colnames(cell_proportions_all_cells_structure_agg_p_stromal)[4] <- "Proportion"

cell_proportions_all_cells_structure_agg <- rbind(cell_proportions_all_cells_structure_agg_p_infiltrated,
                                                  cell_proportions_all_cells_structure_agg_p_internal,
                                                  cell_proportions_all_cells_structure_agg_p_external,
                                                  cell_proportions_all_cells_structure_agg_p_stromal)

library(ggplot2)
cell_proportions_all_cells_structure_agg$Stage <- factor(cell_proportions_all_cells_structure_agg$Stage,
                                                         levels=c("Non-diabetic", "Onset", "Long-duration"))
cell_proportions_all_cells_structure_agg$Location <- factor(cell_proportions_all_cells_structure_agg$Location,
                                                            levels=c("Infiltrated","Internal","External","Stromal"))

cell_proportions_all_cells_structure_agg <- cell_proportions_all_cells_structure_agg[cell_proportions_all_cells_structure_agg$Cell.Type != "All_cells_of_interest",]


cell_proportions_across_images_agg_p_infiltrated <- aggregate(P.Infiltrated.CoI ~ Patient+Stage+Cell.Type, cell_proportions_across_images, mean)
cell_proportions_across_images_agg_p_internal <- aggregate(P.Internal.Margin.CoI ~ Patient+Stage+Cell.Type, cell_proportions_across_images, mean)
cell_proportions_across_images_agg_p_external <- aggregate(P.External.Margin.CoI ~ Patient+Stage+Cell.Type, cell_proportions_across_images, mean)
cell_proportions_across_images_agg_p_stromal <- aggregate(P.Stromal.CoI ~ Patient+Stage+Cell.Type, cell_proportions_across_images, mean)

cell_proportions_across_images_agg_p_infiltrated$Location <- "Infiltrated"
cell_proportions_across_images_agg_p_internal$Location <- "Internal"
cell_proportions_across_images_agg_p_external$Location <- "External"
cell_proportions_across_images_agg_p_stromal$Location <- "Stromal"

colnames(cell_proportions_across_images_agg_p_infiltrated)[4] <- "Proportion"
colnames(cell_proportions_across_images_agg_p_internal)[4] <- "Proportion"
colnames(cell_proportions_across_images_agg_p_external)[4] <- "Proportion"
colnames(cell_proportions_across_images_agg_p_stromal)[4] <- "Proportion"

cell_proportions_across_images_agg <- rbind(cell_proportions_across_images_agg_p_infiltrated,
                                            cell_proportions_across_images_agg_p_internal,
                                            cell_proportions_across_images_agg_p_external,
                                            cell_proportions_across_images_agg_p_stromal)

cell_proportions_across_images_agg$Stage <- factor(cell_proportions_across_images_agg$Stage,
                                                   levels=c("Non-diabetic", "Onset", "Long-duration"))
cell_proportions_across_images_agg$Location <- factor(cell_proportions_across_images_agg$Location,
                                                      levels=c("Infiltrated","Internal","External","Stromal"))

cell_proportions_across_images_agg <- cell_proportions_across_images_agg[cell_proportions_across_images_agg$Cell.Type != "All_cells_of_interest",]


## Calculate p-values for endothelial cells ####
library(clinfun)
trend_p <- vector()
temp <- cell_proportions_across_images_agg[cell_proportions_across_images_agg$Cell.Type == "endothelial",]
for(location in unique(temp$Location)){
  temp2 <- temp[temp$Location == location,]
  g <- factor(temp2$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
  p1 <- jonckheere.test(temp2$Proportion, g, alternative = c("increasing"))$p.value
  p2 <- jonckheere.test(temp2$Proportion, g, alternative = c("decreasing"))$p.value
  trend_p <- rbind(trend_p, c(location, p1, p2))
}
trend_p <- data.frame(trend_p)
colnames(trend_p) <- c("Location", "p_increasing", "p_decreasing")

trend_p$p_increasing <- as.numeric(trend_p$p_increasing)
trend_p$p_decreasing <- as.numeric(trend_p$p_decreasing)

cell_proportions_across_images_agg$Location <- factor(cell_proportions_across_images_agg$Location,
                                                      levels=c("Stromal", "External", "Internal", "Infiltrated"))
ggplot(subset(cell_proportions_across_images_agg, Cell.Type %in% c("endothelial")), 
       aes(x=Stage, y=Proportion))+
  geom_boxplot(aes(fill=Stage))+
  facet_grid(Location ~ ., scale="free_y")+
  ggtitle("Endothelial")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


## Composition of islets ####
load("~/Documents/paper/paper_new_figs/Objects/cell_proportions_in_islet_all.Rdata")
cell_proportions_in_islet_all <- as.data.frame(cell_proportions_in_islet_all)
cell_proportions_in_islet_all$Stage <- donors$stage[match(cell_proportions_in_islet_all$Slide, donors$slide)]

###Looking at the proportion of islet cells in the islets
##Percentage inside vs. the border
####Looking at non-islet cells in the islets
cell_proportions_in_islet_all_only_islet <- cell_proportions_in_islet_all[cell_proportions_in_islet_all$Cell.Type %in% c("alpha",  "beta", "delta", "gamma"),]
cell_proportions_in_islet_all_only_islet$Border <- as.numeric(as.character(cell_proportions_in_islet_all_only_islet$Border))
cell_proportions_in_islet_all_only_islet$Inside <- as.numeric(as.character(cell_proportions_in_islet_all_only_islet$Inside))
cell_proportions_in_islet_all_only_islet$Outside <- as.numeric(as.character(cell_proportions_in_islet_all_only_islet$Outside))

cell_proportions_in_islet_all_only_islet$Inside_plus_border <- cell_proportions_in_islet_all_only_islet$Inside + cell_proportions_in_islet_all_only_islet$Border

cell_proportions_in_islet_all_only_islet$Total_cell_type <- apply(cell_proportions_in_islet_all_only_islet[,c(2:4)], 1, sum) # should be 2:4

total_inside <- aggregate(Inside_plus_border ~ Image, cell_proportions_in_islet_all_only_islet, sum)

cell_proportions_in_islet_all_only_islet$Total_inside <- total_inside$Inside[match(cell_proportions_in_islet_all_only_islet$Image,
                                                                                   total_inside$Image)]

total_image <- aggregate(Total_cell_type ~ Image, cell_proportions_in_islet_all_only_islet, sum)

cell_proportions_in_islet_all_only_islet$Total_image <- total_image$Total_cell_type[match(cell_proportions_in_islet_all_only_islet$Image,
                                                                                          total_image$Image)]

cell_proportions_in_islet_all_only_islet$Percentage_by_cell_type_inside <- (cell_proportions_in_islet_all_only_islet$Inside/cell_proportions_in_islet_all_only_islet$Total_cell_type)*100

cell_proportions_in_islet_all_only_islet$Percentage_by_area_inside <- (cell_proportions_in_islet_all_only_islet$Inside/cell_proportions_in_islet_all_only_islet$Total_inside)*100

cell_proportions_in_islet_all_only_islet$Percentage_by_image_inside <- (cell_proportions_in_islet_all_only_islet$Inside/cell_proportions_in_islet_all_only_islet$Total_image)*100

## By percentage calculate by area ####
cell_proportions_in_islet_all_only_islet_area <- cell_proportions_in_islet_all_only_islet[,c("Patient", "Slide", "Image", "Stage", "Cell.Type",
                                                                                             #"Percentage_by_area_border",
                                                                                             #"Percentage_by_area_outside",
                                                                                             "Percentage_by_area_inside")]
temp_inside <- aggregate(Percentage_by_area_inside ~ Patient+Stage+Cell.Type, cell_proportions_in_islet_all_only_islet_area, mean)
cell_proportions_in_islet_all_only_islet_area2 <- data.frame(temp_inside)

library(reshape2)
cell_proportions_in_islet_all_only_islet_area <- melt(cell_proportions_in_islet_all_only_islet_area2)
cell_proportions_in_islet_all_only_islet_area$variable <- as.character(cell_proportions_in_islet_all_only_islet_area$variable)
cell_proportions_in_islet_all_only_islet_area$variable[cell_proportions_in_islet_all_only_islet_area$variable == "Percentage_by_area_inside"] <- "Inside"
colnames(cell_proportions_in_islet_all_only_islet_area)[4:5] <- c("Area", "Percentage_by_area")
cell_proportions_in_islet_all_only_islet_area$Stage <- factor(cell_proportions_in_islet_all_only_islet_area$Stage,
                                                              levels=c("Non-diabetic","Onset","Long-duration"))
##PAPER 
ggplot(subset(cell_proportions_in_islet_all_only_islet_area, Area == "Inside"), aes(x=Stage, y=Percentage_by_area))+
  geom_boxplot(aes(fill=Stage))+
  facet_grid(.~Cell.Type, scales='free_x')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


library(clinfun)
temp <- subset(cell_proportions_in_islet_all_only_islet_area, Area == "Inside")
temp_alpha <- subset(temp, Cell.Type == "alpha")
temp_beta <- subset(temp, Cell.Type == "beta")
temp_delta <- subset(temp, Cell.Type == "delta")
temp_gamma <- subset(temp, Cell.Type == "gamma")

g <- factor(temp_alpha$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
jonckheere.test(temp_alpha$Percentage_by_area, g, alternative = c("increasing"))

g <- factor(temp_beta$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
jonckheere.test(temp_beta$Percentage_by_area, g, alternative = c("decreasing"))

g <- factor(temp_delta$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
jonckheere.test(temp_delta$Percentage_by_area, g, alternative = c("increasing"))

g <- factor(temp_gamma$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
jonckheere.test(temp_gamma$Percentage_by_area, g, alternative = c("increasing"))


## Islet size ####
##Total cells in an image

temp <- cell_proportions_in_islet_all
temp$Border <- as.numeric(as.character(temp$Border))
temp$Inside <- as.numeric(as.character(temp$Inside))
temp$Outside <- as.numeric(as.character(temp$Outside))
temp$Total_per_cell_type <- apply(temp[,c(2:4)], 1, sum)
total_per_image <- aggregate(Total_per_cell_type ~ Patient+Slide+Image+Stage, temp, sum)

temp_islet_cells <- temp[temp$Cell.Type %in% c("alpha", "beta", "delta", "gamma"),]
temp_islet_cells$Total_in_islets_islet_cells <- apply(temp_islet_cells[,c(2:3)], 1, sum)

total_islets_islet_cells <- aggregate(Total_in_islets_islet_cells ~ Patient+Slide+Image+Stage, temp_islet_cells, sum)

totals <- total_per_image
colnames(totals)[5] <- "Total_in_image"

totals$Total_in_islets_islet_cells <- total_islets_islet_cells$Total_in_islets_islet_cells[match(totals$Image, total_islets_islet_cells$Image)]


totals$Percentage_in_islets_islet_cells <- (totals$Total_in_islets_islet_cells/totals$Total_in_image)*100
totals$Stage <- factor(totals$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"))

totals_islet_cells <- aggregate(Percentage_in_islets_islet_cells ~ Patient + Stage, totals, mean)

totals_islet_cells <- aggregate(Total_in_islets_islet_cells ~ Patient + Stage, totals, mean)

ggplot(totals_islet_cells, aes(x=Stage, y=Total_in_islets_islet_cells))+
  geom_boxplot(aes(fill=Stage))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

library(clinfun)

g <- factor(totals_islet_cells$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
jonckheere.test(totals_islet_cells$Total_in_islets_islet_cells, g, alternative = c("decreasing"))

## Distance to border ####
load("~/Documents/paper/paper_new_figs/Objects/distances_to_border_all.Rdata")
load("~/Documents/paper/paper_new_figs/Objects/proportion_cells_all.Rdata")

##Only consider cases where there are at least 10? cells of the type in the image?

##Find in which images we have at least 10 cells of a type
proportion_cells_all_filtered <- proportion_cells_all[proportion_cells_all$Number_of_celltype >= 10,]

distances_to_border_all$Stage <- donors$stage[match(distances_to_border_all$Slide, donors$slide)]

distances_to_border_all_filtered <- distances_to_border_all[distances_to_border_all$Image %in% proportion_cells_all_filtered$Image,]
##The NAs is when there aren't enough cells of that type in a region. Hence the cutoff of 10 is a bit arbitrary

distances_to_border_all_filtered <- distances_to_border_all_filtered[!is.na(distances_to_border_all_filtered$St.dev_d),]

distances_to_border_agg <- aggregate(Mean_d ~ Cell.Type+Area+Patient+Stage, distances_to_border_all_filtered, mean)

distances_to_border_agg$Stage <- factor(distances_to_border_agg$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"))

distances_to_border_agg <- distances_to_border_agg[distances_to_border_agg$Cell.Type != "All_cell_types_of_interest",]


library(clinfun)
distances_to_border_agg_endothelial <- subset(distances_to_border_agg, Cell.Type=="endothelial")
distances_to_border_agg_stroma <- subset(distances_to_border_agg_endothelial, Area == "Stroma")
distances_to_border_agg_islet <- subset(distances_to_border_agg_endothelial, Area == "Within_border_area")

g <- factor(distances_to_border_agg_islet$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
jonckheere.test(distances_to_border_agg_islet$Mean_d, g, alternative = c("decreasing"))

g <- factor(distances_to_border_agg_stroma$Stage, levels=c("Non-diabetic", "Onset", "Long-duration"), ordered=TRUE)
jonckheere.test(distances_to_border_agg_stroma$Mean_d, g, alternative = c("increasing"))

#PAPER
ggplot(distances_to_border_agg_endothelial, aes(x=Stage, y=Mean_d))+
  geom_boxplot(aes(fill=Stage))+
  ylab("Distance to margin")+
  facet_grid(Area~., scale="free_y")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))




## plot one example image ####
load("data/Objects/real.Rda")

image <- Real[[232]]
categories <- c("alpha","beta","delta","gamma", "ductal","acinar",
                "Th", "otherimmune", "Tc","macrophage", "naiveTc","neutrophil","B", 
                "endothelial","stromal", "unknown")
names <- c(rep("Islet", 4), rep("Exorine",2), rep("Immune", 7), "Endothelial", rep("Others", 2))
image <- define_celltypes(image, categories = categories, category_colname = "Cell.Type",
                          names = names, new_colname = "Cell.Type4")
image$Cell.Type4 <- factor(image$Cell.Type4, levels = c("Islet", "Exorine", "Immune", "Endothelial", "Others"))
new_categories <- levels(colData(image)$Cell.Type4)
colour_vector <- c("brown", "orange", "purple", "black", "gray" )
plot_cell_categories(image, new_categories, colour_vector, "Cell.Type4")

image_border <- identify_bordering_cells(image, "Islet", feature_colname = "Cell.Type4", ahull_alpha = 30)
image_dist <- calculate_distance_to_margin(image_border)
image_structure <- define_structure(image_dist, cell_types_of_interest = "Endothelial", feature_colname = "Cell.Type4")
plot_cell_categories(image_structure, feature_colname = "Structure")

max(spatialCoords(image)[,1]) #685
max(spatialCoords(image)[,2]) #540.625

## compute number of islets ####
library(SPIAT)
library(dplyr)
library(tibble)

# FUNCTIONS ----
get_colData <- function(spe_object){
  formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
  formatted_data <- cbind(formatted_data, 
                          data.frame(SpatialExperiment::spatialCoords(spe_object)))
  if (is.null(formatted_data$Cell.ID)){
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
  }
  
  # delete column `sample_id`
  formatted_data$sample_id <- NULL
  
  return(formatted_data)
}
min_min_dist <- function(spe_object) {
  formatted_data <- get_colData(spe_object)
  
  #extract the cell coordinates
  all_cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]
  
  #CHECK
  if (nrow(all_cell_cords) == 0) {
    stop("No cells found in average minimum distance calculation")
  }
  
  #calculate the closest 2 neighbours, 1st being itself
  all_closest <- RANN::nn2(data = all_cell_cords, k = 2)
  
  #grab the distances and find the average
  all_closest_dist <- all_closest$nn.dists[,2]
  min_min_distance <- min(all_closest_dist)
  
  return(min_min_distance)
}
preprocess <- function(){
  items <- list()
  load("data.Rdata")
  # 12 patients, 24 slides, 845 cores, 
  metadata <- read.delim("Version_1/Metadata.csv",sep=",")
  panel <- read.delim("Version_1/Panel.csv",sep=",")
  celltypes <- read.delim("Version_2/CellTypes.csv",sep=",")
  donors <- read.delim("Version_2/Donors.csv",sep=",")
  
  celltypes$cell <- substr(celltypes$id, 5, nchar(celltypes$id))
  celltypes$region <- substr(celltypes$id, 1, 1)
  celltypes$patient <- donors$case[match(celltypes$region, donors$slide)]
  
  image_map <- unique(celltypes[,c("core", "region", "patient")])
  image_map$image <- 1:nrow(image_map)
  image_map$Image_ID <- paste(image_map$patient,
                              image_map$region,
                              image_map$image, sep="_")
  
  celltypes$Image_ID <- image_map$Image_ID[match(celltypes$core, image_map$core)]
  celltypes$CellType2 <- celltypes$CellType
  celltypes$CellType2[celltypes$CellType2 %in% c("alpha", "delta", "beta",  "gamma")] <- "islet"
  ## QUESTION: How do you know the IDs match?
  data$Image_ID <- image_map$Image_ID[match(data$ImageNumber, image_map$image)]
  
  # return
  for (item in c("data", "celltypes")){
    items[[item]] <- eval(parse(text = item))
  }
  
  return(items)
}
### get the ratio of the image sides
get_window_ratio <- function(args, metrics){
  args[, "window_ratio"] <- as.numeric(metrics[, "Y_max"])/as.numeric(metrics[, "X_max"])
  return(args)
}
### calculate the proportions
# this function need to return a list with one arg table and two index tables
calculate_proportions <- function(args, metrics){
  args[, "n_infil_beta"] <- as.numeric(metrics[, "beta_Inside"]) + as.numeric(metrics[, "beta_Border"]) 
  args[, "n_infil_nonbeta"] <- as.numeric(metrics[, "alpha_Inside"]) + 
    as.numeric(metrics[, "gamma_Inside"]) + as.numeric(metrics[, "delta_Inside"]) +
    as.numeric(metrics[, "alpha_Border"]) + 
    as.numeric(metrics[, "gamma_Border"]) + as.numeric(metrics[, "delta_Border"])
  
  args[, "n_infil_immune"] <- 
    as.numeric(metrics[, "Th_Inside"]) + 
    as.numeric(metrics[, "Tc_Inside"]) + 
    as.numeric(metrics[, "macrophage_Inside"]) +
    as.numeric(metrics[, "neutrophil_Inside"]) + 
    as.numeric(metrics[, "otherimmune_Inside"]) + 
    as.numeric(metrics[, "B_Inside"]) + 
    as.numeric(metrics[, "naiveTc_Inside"]) +
    as.numeric(metrics[, "Th_Border"]) + 
    as.numeric(metrics[, "Tc_Border"]) + as.numeric(metrics[, "macrophage_Border"]) +
    as.numeric(metrics[, "neutrophil_Border"]) + as.numeric(metrics[, "otherimmune_Border"]) + 
    as.numeric(metrics[, "B_Border"]) + as.numeric(metrics[, "naiveTc_Border"])
  
  args[, "n_infil_endothelial"] <- as.numeric(metrics[, "endothelial_Inside"]) +
    as.numeric(metrics[, "endothelial_Border"]) 
  args[, "n_infil_others"] <- as.numeric(metrics[, "acinar_Inside"]) + 
    as.numeric(metrics[, "ductal_Inside"]) + as.numeric(metrics[, "unknown_Inside"]) + 
    as.numeric(metrics[, "stromal_Inside"]) +as.numeric(metrics[, "acinar_Border"]) + 
    as.numeric(metrics[, "ductal_Border"]) + as.numeric(metrics[, "unknown_Border"]) + 
    as.numeric(metrics[, "stromal_Border"])
  
  args[, "n_infil_all"] <- args[, "n_infil_beta"] + args[, "n_infil_nonbeta"] +
    args[, "n_infil_immune"] + args[, "n_infil_endothelial"] + args[, "n_infil_others"]
  
  args[, "n_stromal_beta"] <- as.numeric(metrics[, "beta_Outside"])
  args[, "n_stromal_nonbeta"] <- as.numeric(metrics[, "alpha_Outside"]) + 
    as.numeric(metrics[, "gamma_Outside"]) + as.numeric(metrics[, "delta_Outside"])
  args[, "n_stromal_immune"] <- 
    as.numeric(metrics[, "Th_Outside"]) + 
    as.numeric(metrics[, "Tc_Outside"]) + 
    as.numeric(metrics[, "macrophage_Outside"]) +
    as.numeric(metrics[, "neutrophil_Outside"]) + 
    as.numeric(metrics[, "otherimmune_Outside"]) + 
    as.numeric(metrics[, "B_Outside"]) + 
    as.numeric(metrics[, "naiveTc_Outside"]) 
  args[, "n_stromal_endothelial"] <- as.numeric(metrics[, "endothelial_Outside"]) 
  args[, "n_stromal_others"] <- as.numeric(metrics[, "acinar_Outside"]) + 
    as.numeric(metrics[, "ductal_Outside"]) + as.numeric(metrics[, "unknown_Outside"]) + 
    as.numeric(metrics[, "stromal_Outside"])
  args[, "n_stromal_all"] <- args[, "n_stromal_beta"] + args[, "n_stromal_nonbeta"] +
    args[, "n_stromal_immune"] + args[, "n_stromal_endothelial"] + args[, "n_stromal_others"]
  
  for (type in c("beta", "nonbeta", "immune", "endothelial", "others")){
    args[, paste("p_infil", type, sep = "_")] <-  
      round(args[, paste("n_infil", type, sep = "_")]/args[, "n_infil_all"], 2)
    args[, paste("p_stromal", type, sep = "_")] <-  
      round(args[, paste("n_stromal", type, sep = "_")]/args[, "n_stromal_all"], 2)
  }
  args[, "p_infil_idx"] <- seq_len(nrow(args))
  args[, "p_stromal_idx"] <- seq_len(nrow(args))
  infil_idx_t <- args[, c("p_infil_idx", "p_infil_beta", "p_infil_nonbeta", 
                          "p_infil_immune", "p_infil_endothelial", "p_infil_others")]
  stromal_idx_t <- args[, c("p_stromal_idx", "p_stromal_beta", "p_stromal_nonbeta", 
                            "p_stromal_immune", "p_stromal_endothelial", "p_stromal_others")]
  args[, "n_cluster_cells"] <- args[, "n_infil_all"]
  args[, c("p_infil_beta", "p_infil_nonbeta", 
           "p_infil_immune", "p_infil_endothelial", "p_infil_others",
           "n_infil_beta", "n_infil_nonbeta", "n_infil_immune", 
           "n_infil_endothelial", "n_infil_others", "n_infil_all",
           "p_stromal_beta", "p_stromal_nonbeta", 
           "p_stromal_immune", "p_stromal_endothelial", "p_stromal_others",
           "n_stromal_beta", "n_stromal_nonbeta", "n_stromal_immune", 
           "n_stromal_endothelial", "n_stromal_others", "n_stromal_all")] <- NULL
  
  return(list(args = args, infil_idx_t = infil_idx_t, stromal_idx_t = stromal_idx_t))
}
### separate the min_min_dist column (need to shuffle together with avg_min_dist)
separate_min_min_dist <- function(args_list){
  args <- args_list$args
  args[, "dist_idx"] <- seq_len(nrow(args))
  dist_idx_t <- args[, c("dist_idx", "avg_min_dist", "min_min_dist")]
  args[, c("avg_min_dist", "min_min_dist")] <- NULL
  args_list$args <- args
  args_list$dist_idx_t <- dist_idx_t
  return(args_list)
}

# PREPROCESS ####
items <- preprocess()
data <- items$data
celltypes <- items$celltypes
rm(items)

# LOOP get metrics ####
n_cells <- image_size <- proportion_cells <- overall_min_distance <- 
  min_min_distance <- n_clusters <- cell_proportions_in_islet <- list()
pb <- txtProgressBar(style = 3)
for(j in seq_len(length(unique(data$Image_ID)))){
  # for (j in c(415, 522, 566, 654, 687, 688, 693, 781)){ #some bordering cell detection do not return bordering cells
  i <- unique(data$Image_ID)[j]
  pdf(paste("Plots/Plots_", i, ".pdf", sep=""), width=8)
  temp_image <- data[data$Image_ID == i,]
  temp_celltypes <- celltypes[celltypes$Image_ID == i,]
  
  temp_image$CellType <- temp_celltypes$CellType[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  temp_image$CellType2 <- temp_celltypes$CellType2[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  temp_image$CellCat <- temp_celltypes$CellCat[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  
  g <- ggplot(temp_image, aes(x=Location_Center_X, y=Location_Center_Y))+
    geom_point(aes(colour=CellType))
  print(g)
  g <- ggplot(temp_image, aes(x=Location_Center_X, y=Location_Center_Y))+
    geom_point(aes(colour=CellCat))
  print(g)
  
  map <- data.frame(CellType = temp_image$CellType,
                    CellType2 = temp_image$CellType2,
                    CellCat = temp_image$CellCat)
  map <- unique(map)
  
  formatted_image <- 
    format_image_to_spe(format = "general", phenotypes = temp_image$CellType, 
                        coord_x = temp_image$Location_Center_X, coord_y = temp_image$Location_Center_Y)
  rownames(colData(formatted_image)) <- paste("Cell", 1:nrow(colData(formatted_image)), sep="_")
  
  formatted_image <- 
    define_celltypes(formatted_image,  categories = map$CellType, category_colname = "Phenotype", 
                     names = map$CellCat,new_colname = "Cell.Cat")
  
  formatted_image <- 
    define_celltypes(formatted_image, categories = map$CellType, category_colname = "Phenotype", 
                     names = map$CellType,new_colname = "Cell.Type")
  
  formatted_image <- 
    define_celltypes( formatted_image,  categories = map$CellType,  category_colname = "Phenotype", 
                      names = map$CellType2, new_colname = "Cell.Type2")
  
  #Total number of cells
  n_cells[[i]] <- ncol(formatted_image)
  
  #Image size
  image_size[[i]] <- c(X_coord = max(spatialCoords(formatted_image)[,1]), 
                       Y_coord = max(spatialCoords(formatted_image)[,2]))
  
  ##Proportion of cells in the entire image. You can get islet size from here (see my code for the paper below)
  proportion_cells[[i]] <- calculate_cell_proportions(formatted_image, 
                                                      reference_celltypes=NULL, 
                                                      feature_colname ="Cell.Type",
                                                      plot.image = FALSE)
  # average_min_dist
  overall_min_distance[[i]] <- average_minimum_distance(formatted_image)
  # min_min_dist
  min_min_distance[[i]] <- min_min_dist(formatted_image)
  
  cells <- unique(colData(formatted_image)$Cell.Type)
  cells <- cells[!(cells %in% c("unknown", "otherimmune"))]
  cells <- c(cells, "islet")
  
  ##Identify margin cells####
  formatted_border <- tryCatch(
    identify_bordering_cells(formatted_image, reference_cell = "islet", 
                             ahull_alpha = 20, feature_colname="Cell.Type2", 
                             n_of_polygons = 1, plot_final_border = TRUE),
    error=function(e) e, warning=function(w) w)
  
  if (is(formatted_border, "warning")){
    formatted_border <- tryCatch(
      identify_bordering_cells(formatted_image, reference_cell = "islet", 
                               ahull_alpha = 30, feature_colname="Cell.Type2", 
                               n_of_polygons = 1, plot_final_border = TRUE),
      error=function(e) e, warning=function(w) w)
  }
  
  if (is(formatted_border, "warning")){
    formatted_border <- tryCatch(
      identify_bordering_cells(formatted_image, reference_cell = "islet", 
                               ahull_alpha = 40, feature_colname="Cell.Type2", 
                               n_of_polygons = 1, plot_final_border = TRUE),
      error=function(e) e, warning=function(w) w)
  }
  
  if (is(formatted_border, "warning")){
    formatted_border <- tryCatch(
      identify_bordering_cells(formatted_image, reference_cell = "islet", 
                               ahull_alpha = 40, feature_colname="Cell.Type2", 
                               n_of_polygons = 1, plot_final_border = TRUE, 
                               n_to_exclude = 5),
      error=function(e) e, warning=function(w) w)
  }
  
  
  n_clusters[[i]] <- attr(formatted_border, "n_of_clusters")
  
  formatted_distance <- calculate_distance_to_margin(formatted_border)
  
  ##Cells in structure
  non_islet <- c("acinar","ductal","stromal","endothelial","Th","macrophage","neutrophil",
                 "Tc","naiveTc","B")
  islet <- c("alpha", "delta", "beta", "gamma")
  
  ##Composition of islets
  formatted_structure <- define_structure(
    formatted_distance, cell_types_of_interest = islet, 
    feature_colname = "Cell.Type", n_margin_layers = 0)
  plot_cell_categories(formatted_structure, feature_colname = "Structure")
  
  cell_proportions_in_islet[[i]] <- table(formatted_structure@colData[,c("Cell.Type", "Region")])
  
  # save objects
  save(n_cells, file=paste("Objects/n_cells.Rdata", sep=""))  
  save(image_size, file=paste("Objects/image_size.Rdata", sep=""))  
  save(proportion_cells, file=paste("Objects/proportion_cells.Rdata", sep=""))
  save(overall_min_distance, file=paste("Objects/overall_min_distance.Rdata", sep=""))
  save(overall_min_distance, file=paste("Objects/min_min_distance.Rdata", sep=""))
  save(n_clusters, file=paste("Objects/n_clusters.Rdata", sep=""))
  save(cell_proportions_in_islet, file=paste("Objects/cell_proportions_in_islet.Rdata", sep=""))
  print(i)
  dev.off()
  setTxtProgressBar(pb, j/length(unique(data$Image_ID)))
}
close(pb)
rm(g); rm(temp_image); rm(temp_celltypes); rm(formatted_structure)
rm(formatted_border); rm(formatted_image); rm(formatted_distance)
rm(map); rm(i); rm(j)

# Collapsing metrics to single data frames ####
load(paste("Objects/n_cells.Rdata", sep=""))
load(paste("Objects/image_size.Rdata", sep=""))  
load(paste("Objects/overall_min_distance.Rdata", sep=""))
load(paste("Objects/min_min_distance.Rdata", sep=""))
load(paste("Objects/proportion_cell.Rdata", sep=""))
load(paste("Objects/cell_proportions_in_islet.Rdata", sep=""))

agg_metrics <- data.frame(matrix(ncol = 25, nrow = 0))

for(image_ID in unique(data$Image_ID)){
  patient <- strsplit(image_ID, "_")[[1]][1]
  slide <- strsplit(image_ID, "_")[[1]][2]
  image <- strsplit(image_ID, "_")[[1]][3]
  
  # image id, n_cells, image size, min dist, n_clusters
  temp_row <- c(Patient = patient, Slide = slide, Image = image, 
                n_cells = as.numeric(as.character(n_cells[[image_ID]])),
                X_max = ceiling(as.numeric(as.character(image_size[[image_ID]][1]))),
                Y_max = ceiling(as.numeric(as.character(image_size[[image_ID]][2]))),
                avg_min_dist = round(as.numeric(as.character(overall_min_distance[[image_ID]])), 1),
                min_min_dist = round(as.numeric(as.character(min_min_distance[[image_ID]])), 1),
                n_clusters = as.numeric(as.character(n_clusters[[image_ID]])))
  
  # cell count in whole image
  temp <- proportion_cells[[image_ID]]
  temp <- remove_rownames(temp)
  temp <- column_to_rownames(temp, "Cell_type")
  temp[,c("Proportion", "Proportion_name", "Percentage")] <- NULL
  temp <- as.data.frame(t(temp))
  temp_proportion_row <- c(rep(NA, 16))
  names(temp_proportion_row) <- unique(celltypes$CellType)
  temp_proportion_row[names(temp)] <- temp
  temp_row <- unlist(c(temp_row, temp_proportion_row))
  
  # cell count in, on and out of islet border
  for (region in c("Inside", "Border", "Outside")){
    temp <- cell_proportions_in_islet[[image_ID]][, region]
    temp_proportion_islet_row <- c(rep(NA, 16))
    names(temp_proportion_islet_row) <- unique(celltypes$CellType)
    temp_proportion_islet_row[names(temp)] <- temp
    names(temp_proportion_islet_row) <- paste(names(temp_proportion_islet_row), region, sep = "_")
    temp_row <- unlist(c(temp_row, temp_proportion_islet_row))
  }
  
  # aggregate
  agg_metrics <- as.data.frame(rbind(agg_metrics, temp_row))
  
  save(agg_metrics, file = "Objects/agg_metrics.Rdata")
  
  print(image)
}
x <- c("Patient", "Slide", "Image", "n_cells", "X_max",
       "Y_max", "avg_min_dist", "min_min_dist", "n_clusters", unique(celltypes$CellType),
       paste(unique(celltypes$CellType), "Inside", sep = "_"),
       paste(unique(celltypes$CellType), "Border", sep = "_"),
       paste(unique(celltypes$CellType), "Outside", sep = "_"))
colnames(agg_metrics) <- x

agg_metrics[is.na(agg_metrics)] <- 0
save(agg_metrics, file = "Objects/agg_metrics.Rdata")

load("Objects/agg_metrics.Rdata")
df_nclusters <- agg_metrics[, c("Image", "n_clusters")]

metadata <- read.delim("data/Version_1/Metadata.csv",sep=",")
metadata$Image <- 1:845

df_nclusters <- merge(df_nclusters, metadata, by = "Image")

df_nclusters$n_clusters <- as.numeric(df_nclusters$n_clusters)

sum(df_nclusters$n_clusters)

sum(df_nclusters[which(df_nclusters$stage == "Non-diabetic"), "n_clusters"])/
  table(df_nclusters$stage)["Non-diabetic"]
sum(df_nclusters[which(df_nclusters$stage == "Onset"), "n_clusters"])/
  table(df_nclusters$stage)["Onset"]
sum(df_nclusters[which(df_nclusters$stage == "Long-duration"), "n_clusters"])/
  table(df_nclusters$stage)["Long-duration"]
