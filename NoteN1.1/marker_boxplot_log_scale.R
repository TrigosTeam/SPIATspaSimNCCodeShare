bind_info <- function(spe_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
    formatted_data <- cbind(formatted_data,
                            data.frame(SpatialExperiment::spatialCoords(spe_object)))
    #convert rowname to column
    if (is.null(formatted_data$Cell.ID)){
        formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
    }
    # get the intensity matrix
    intensity_matrix <- SummarizedExperiment::assay(spe_object)
    markers <- rownames(intensity_matrix)
    cell_ids <- colnames(intensity_matrix)
    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_t <- data.frame(t(intensity_matrix))
    colnames(intensity_t) <- markers
    
    # bind
    formatted_data <- cbind(formatted_data, intensity_t)
    
    # delete column `sample_id`
    formatted_data$sample_id <- NULL
    
    return(formatted_data)
}

marker_boxplot_log_scale <- function (sce_object, marker)
{
    intensity <- NULL
    formatted_data <- bind_info(sce_object)
    intensity_true <- formatted_data[grepl(marker, formatted_data$Phenotype), 
    ]
    intensity_false <- formatted_data[!grepl(marker, formatted_data$Phenotype), 
    ]
    if (nrow(intensity_true) != 0) {
        intensity_true$intensity <- "P"
    }
    else {
        stop(sprintf("There are no cells positive for %s", marker))
    }
    if (nrow(intensity_false) != 0) {
        intensity_false$intensity <- "N"
    }
    else {
        stop(sprintf("There are no cells negative for %s", marker))
    }
    intensity_by_marker <- rbind(intensity_true, intensity_false)
    intensity_by_marker$intensity <- factor(intensity_by_marker$intensity, 
                                            levels = c("P", "N"))
    title <- paste("Level of ", marker, sep = "")
    give.n <- function(x) {
        return(c(y = max(x) + 1, label = length(x)))
    }
    p <- ggplot(intensity_by_marker, aes(x = intensity, 
                                         y = log(intensity_by_marker[, marker])))
    p <- p + geom_boxplot()
    p <- p + labs(title = title, x = "Marker status (Positive/Negative cell)", 
                  y = "Marker level")
    p <- p + stat_summary(fun.data = give.n, geom = "text")
    p <- p + theme_bw()
    print(p)
}
