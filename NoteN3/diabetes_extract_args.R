# This script is for processing the diabetes images and extract arguments for simulations
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
  load("data/data.Rdata")
  # 12 patients, 24 slides, 845 cores, 
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
  temp_image <- data[data$Image_ID == i,]
  temp_celltypes <- celltypes[celltypes$Image_ID == i,]
  
  temp_image$CellType <- temp_celltypes$CellType[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  temp_image$CellType2 <- temp_celltypes$CellType2[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  temp_image$CellCat <- temp_celltypes$CellCat[match(temp_image$ObjectNumber, temp_celltypes$cell)]
  
  
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

# remove redundant variables ####
rm(cell_proportions_in_islet); rm(image_size); rm(min_min_distance)
rm(overall_min_distance); rm(n_cells); rm(n_clusters); rm(proportion_cells)
rm(temp_proportion_row); rm(image); rm(image_ID); rm(patient); rm(region)
rm(slide); rm(temp); rm(temp_proportion_islet_row); rm(temp_row); rm(x)

# Calculate the arguments for simulation ####
load("Objects/agg_metrics.Rdata")

args_sim <- agg_metrics[,c("Patient", "Slide", "Image", "n_cells", 
                           "avg_min_dist", "min_min_dist", "n_clusters")]
args_sim <- get_window_ratio(args_sim, agg_metrics)
args_sim_list <- calculate_proportions(args_sim, agg_metrics)
# separate min_min_dist from avg_min_dist and put it into another index table
args_sim_list <- separate_min_min_dist(args_sim_list)

save(args_sim_list, file = "Objects/independent_sim_args.Rda")
