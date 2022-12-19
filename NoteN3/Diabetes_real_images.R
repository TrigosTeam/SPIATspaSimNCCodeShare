# this script is for reading and formatting the diabetes dataset 
# (including defining the cell types and identifying islet margin).

# Publically available data downloaded from https://data.mendeley.com/datasets/cydmwsfztj/2

#data <- read.delim("/home/atrigos/Spatial_paper/Public_datasets/Diabetes_mass_cyto/All_Cells.csv",sep=",")
#save(data, file="/home/atrigos/Spatial_paper/Public_datasets/Diabetes_mass_cyto/data.Rdata")

# The following file paths are my local file paths (Files not uploaded to GitHub)
load("data/data.Rdata")
metadata <- read.delim("data/Version_1/Metadata.csv",sep=",")
panel <- read.delim("data/Version_1/Panel.csv",sep=",")
celltypes <- read.delim("data/Version_2/CellTypes.csv",sep=",")
donors <- read.delim("data/Version_2/Donors.csv",sep=",")

# Define celltypes
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
celltypes$CellType3 <- celltypes$CellType
celltypes$CellType4 <- celltypes$CellType
celltypes$CellType2[celltypes$CellType2 %in% c("alpha", "delta", "beta",  "gamma")] <- "islet"
celltypes$CellType3[celltypes$CellType3 %in% c("alpha", "delta", "gamma")] <- "nonbeta"
celltypes$CellType4[celltypes$CellType %in% c("alpha", "delta", "beta",  "gamma")] <- "islet"
data$Image_ID <- image_map$Image_ID[match(data$ImageNumber, image_map$image)]

# Format each image into spe object and save the object list
library(SPIAT)
Real <- list()
for(i in unique(data$Image_ID)){
    
    temp_image <- data[data$Image_ID == i,]
    temp_celltypes <- celltypes[celltypes$Image_ID == i,]
    
    temp_image$CellType <- temp_celltypes$CellType[match(temp_image$ObjectNumber, temp_celltypes$cell)]
    temp_image$CellType2 <- temp_celltypes$CellType2[match(temp_image$ObjectNumber, temp_celltypes$cell)]
    temp_image$CellType3 <- temp_celltypes$CellType3[match(temp_image$ObjectNumber, temp_celltypes$cell)]
    temp_image$CellType4 <- temp_celltypes$CellType4[match(temp_image$ObjectNumber, temp_celltypes$cell)]
    temp_image$CellCat <- temp_celltypes$CellCat[match(temp_image$ObjectNumber, temp_celltypes$cell)]
    
    map <- data.frame(CellType = temp_image$CellType,
                      CellType2 = temp_image$CellType2,
                      CellType3 = temp_image$CellType3,
                      CellType4 = temp_image$CellType4,
                      CellCat = temp_image$CellCat)
    map <- unique(map)
    
    formatted_image <- format_image_to_spe(format = "general", 
                                           phenotypes = temp_image$CellType, 
                                           coord_x = temp_image$Location_Center_X,
                                           coord_y = temp_image$Location_Center_Y)
    rownames(colData(formatted_image)) <- paste("Cell", 1:nrow(colData(formatted_image)), sep="_")
    
    for(type in c("Cell.Cat", "Cell.Type", "Cell.Type2", "Cell.Type3", "Cell.Type4")){
        formatted_image <- define_celltypes(
            formatted_image, categories = map$CellType, 
            category_colname = "Phenotype", 
            names = map[[gsub(".", "", type, fixed = TRUE)]],
            new_colname = type)
    }
    
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
    
    attr(formatted_border, "name") <- i
    Real[[i]] <- formatted_border
}
save(Real, file = "Objects/diabetes_real_spe.Rda") 
