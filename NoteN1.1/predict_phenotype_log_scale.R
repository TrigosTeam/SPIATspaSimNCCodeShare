library(dplyr)
library(ggplot2)
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

predict_phenotype_log_scale <- function (sce_object, thresholds = NULL, tumour_marker, baseline_markers, 
          nuclear_marker = NULL, reference_phenotypes = FALSE, markers_to_phenotype = NULL, 
          plot_distribution = TRUE) 
{
    Marker_level <- NULL
    formatted_data <- get_colData(sce_object)
    intensity_matrix <- SummarizedExperiment::assay(sce_object)
    if (is.null(markers_to_phenotype)) {
        markers <- rownames(intensity_matrix)
    }
    else {
        markers <- markers_to_phenotype
        if (sum(markers %in% rownames(intensity_matrix)) != length(markers)) {
            missing <- markers[!(markers %in% rownames(intensity_matrix))]
            stop(sprintf("There are no intensity values for %s in your input dataset. Please check markers\n", 
                         missing))
        }
        intensity_matrix <- intensity_matrix[markers, ]
    }
    if (is.element(tumour_marker, markers) == FALSE) {
        stop("Tumour marker not found")
    }
    cell_ids <- colnames(intensity_matrix)
    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_matrix_t <- t(intensity_matrix)
    intensity_df <- data.frame(intensity_matrix_t)
    colnames(intensity_df) <- markers
    formatted_data <- cbind(formatted_data, intensity_df)
    markers_no_tumour <- markers[markers != tumour_marker]
    if (!is.null(nuclear_marker)) {
        markers_no_tumour <- markers_no_tumour[markers_no_tumour != 
                                                   nuclear_marker]
    }
    selected_valley_xcord <- list()
    for (marker in markers) {
        if ((!is.null(nuclear_marker)) && (marker == nuclear_marker)) {
            actual_phenotype <- data.frame(rep(1, nrow(formatted_data)))
            colnames(actual_phenotype) <- paste(marker, "_actual_phenotype", 
                                                sep = "")
            formatted_data <- cbind(formatted_data, actual_phenotype)
        }
        else {
            cell_IDs <- data.frame(formatted_data$Cell.ID)
            colnames(cell_IDs) <- "Cell.ID"
            actual_phenotype <- data.frame(rep(0, nrow(formatted_data)))
            colnames(actual_phenotype) <- "Phenotype_status"
            actual_phenotype <- cbind(actual_phenotype, cell_IDs)
            rows <- formatted_data[grepl(marker, formatted_data$Phenotype), 
            ]
            if (nrow(rows) > 0) {
                actual_phenotype[actual_phenotype$Cell.ID %in% 
                                     rows$Cell.ID, ]$Phenotype_status <- 1
                actual_phenotype <- data.frame(actual_phenotype[, 
                                                                1])
                colnames(actual_phenotype) <- paste(marker, "_actual_phenotype", 
                                                    sep = "")
                formatted_data <- cbind(formatted_data, actual_phenotype)
            }
        }
    }
    for (marker in markers_no_tumour) {
        marker_specific_level <- formatted_data[, marker]
        if (!is.null(thresholds) && !is.na(thresholds[match(marker, 
                                                            markers)])) {
            marker_threshold <- thresholds[match(marker, markers)]
            sprintf("(%s has threshold specified: %s)", marker, 
                    as.character(marker_threshold))
            selected_valley_xcord[[marker]] <- NULL
            predictions_by_threshold <- data.frame(mmand::threshold(marker_specific_level, 
                                                                    level = marker_threshold))
        }
        else {
            intensity_density <- stats::density(marker_specific_level, 
                                                na.rm = TRUE)
            valleys <- pracma::findpeaks(-(intensity_density)$y)
            valley_ycords <- valleys[, 1] * -1
            index <- match(valley_ycords, intensity_density$y)
            valley_xcords <- intensity_density$x[index]
            valley_df <- data.frame(cbind(valley_xcords, valley_ycords))
            ycord_max_density <- max(intensity_density$y)
            xcord_max_density_index <- match(ycord_max_density, 
                                             intensity_density$y)
            xcord_max_density <- intensity_density$x[xcord_max_density_index]
            density_threshold_for_valley <- 0.25 * ycord_max_density
            valley_df <- valley_df[valley_df$valley_xcords >= 
                                       xcord_max_density, ]
            valley_df <- valley_df[valley_df$valley_ycords <= 
                                       density_threshold_for_valley, ]
            selected_valley_xcord[[marker]] <- valley_df$valley_xcords[1]
        }
        predictions_by_threshold <- data.frame(mmand::threshold(marker_specific_level, 
                                                                level = selected_valley_xcord[[marker]]))
        colnames(predictions_by_threshold) <- paste(marker, "_predicted_phenotype", 
                                                    sep = "")
        formatted_data <- cbind(formatted_data, predictions_by_threshold)
    }
    baseline_cells <- vector()
    for (marker in baseline_markers) {
        temp <- formatted_data[, colnames(formatted_data) %in% 
                                   c("Cell.ID", paste0(marker, "_predicted_phenotype"))]
        temp <- temp[temp[, 2] == 1, ]
        baseline_cells <- c(baseline_cells, temp$Cell.ID)
    }
    baseline_cells <- unique(baseline_cells)
    formatted_data_baseline <- formatted_data[formatted_data$Cell.ID %in% 
                                                  baseline_cells, tumour_marker]
    if (length(unique(formatted_data_baseline)) != 2) {
        cutoff_for_tumour <- stats::quantile(formatted_data_baseline, 
                                             0.95)
    }
    else {
        cutoff_for_tumour <- max(formatted_data_baseline) - min(formatted_data_baseline)/2
    }
    tumour_specific_level <- formatted_data[, tumour_marker]
    if (!is.null(thresholds)) {
        marker_threshold <- thresholds[match(tumour_marker, markers)]
        sprintf("(%s has threshold specified: %s)", tumour_marker, 
                as.character(marker_threshold))
        selected_valley_xcord[[marker]] <- NULL
        predictions_by_threshold <- data.frame(mmand::threshold(tumour_specific_level, 
                                                                level = marker_threshold))
    }
    else {
        intensity_density <- stats::density(tumour_specific_level, 
                                            na.rm = TRUE)
        valleys <- pracma::findpeaks(-(intensity_density)$y)
        valley_ycords <- valleys[, 1] * -1
        index <- match(valley_ycords, intensity_density$y)
        valley_xcords <- intensity_density$x[index]
        valley_df <- data.frame(cbind(valley_xcords, valley_ycords))
        selected_valley_xcord[[tumour_marker]] <- valley_df$valley_xcords[1]
        final_threshold <- ifelse(selected_valley_xcord[[tumour_marker]] < 
                                      cutoff_for_tumour, selected_valley_xcord[[tumour_marker]], 
                                  cutoff_for_tumour)
        selected_valley_xcord[[tumour_marker]] <- final_threshold
        predictions_by_threshold <- data.frame(mmand::threshold(tumour_specific_level, 
                                                                level = final_threshold))
        colnames(predictions_by_threshold) <- paste(tumour_marker, 
                                                    "_predicted_phenotype", sep = "")
        formatted_data <- cbind(formatted_data, predictions_by_threshold)
    }
    predicted_data <- formatted_data
    if (!is.null(nuclear_marker)) {
        markers <- markers[markers != nuclear_marker]
    }
    if (reference_phenotypes) {
        for (marker in markers) {
            print(marker)
            marker_status_name <- paste(marker, "_status", sep = "")
            marker_actual_exp_colname <- paste(marker, "_actual_phenotype", 
                                               sep = "")
            marker_pred_exp_colname <- paste(marker, "_predicted_phenotype", 
                                             sep = "")
            accuracy_df <- data.frame(rep(NA, nrow(predicted_data)))
            colnames(accuracy_df) <- "status"
            if (marker_actual_exp_colname %in% colnames(predicted_data)) {
                marker_exp_actual_pred <- predicted_data[, c(marker_actual_exp_colname, 
                                                             marker_pred_exp_colname)]
                colnames(marker_exp_actual_pred) <- c("actual", 
                                                      "pred")
                accuracy_df <- cbind(accuracy_df, marker_exp_actual_pred)
                if (nrow(accuracy_df[accuracy_df$actual == 1 & 
                                     accuracy_df$pred == 1, ])) {
                    accuracy_df[accuracy_df$actual == 1 & accuracy_df$pred == 
                                    1, ]$status <- "TP"
                }
                if (nrow(accuracy_df[accuracy_df$actual == 0 & 
                                     accuracy_df$pred == 0, ])) {
                    accuracy_df[accuracy_df$actual == 0 & accuracy_df$pred == 
                                    0, ]$status <- "TN"
                }
                if (nrow(accuracy_df[accuracy_df$actual == 0 & 
                                     accuracy_df$pred == 1, ])) {
                    accuracy_df[accuracy_df$actual == 0 & accuracy_df$pred == 
                                    1, ]$status <- "FP"
                }
                if (nrow(accuracy_df[accuracy_df$actual == 1 & 
                                     accuracy_df$pred == 0, ])) {
                    accuracy_df[accuracy_df$actual == 1 & accuracy_df$pred == 
                                    0, ]$status <- "FN"
                }
                accuracy_df <- data.frame(accuracy_df[, "status"])
                colnames(accuracy_df) <- "status"
                marker_specific_level <- predicted_data[, marker]
                marker_specific_level <- data.frame(marker_specific_level)
                colnames(marker_specific_level) <- "Marker_level"
                level_and_accuracy <- cbind(marker_specific_level, 
                                            accuracy_df)
                TP_count <- nrow(level_and_accuracy[level_and_accuracy$status == 
                                                        "TP", ])
                TN_count <- nrow(level_and_accuracy[level_and_accuracy$status == 
                                                        "TN", ])
                FP_count <- nrow(level_and_accuracy[level_and_accuracy$status == 
                                                        "FP", ])
                FN_count <- nrow(level_and_accuracy[level_and_accuracy$status == 
                                                        "FN", ])
                sprintf("For %s:", marker)
                sprintf("TP:%i TN:%i FP:%i FN:%i", TP_count, 
                        TN_count, FP_count, FN_count)
                if (plot_distribution) {
                    p <- ggplot(level_and_accuracy, aes(x = log(Marker_level))) + 
                        geom_density()
                    title <- paste("Density distribution of", marker, 
                                   sep = " ")
                    p <- p + labs(title = title, x = "Level of intensity (log)", 
                                  y = "Density")
                    if (!is.null(selected_valley_xcord[[marker]])) {
                        p <- p + geom_vline(aes(xintercept = log(selected_valley_xcord[[marker]])), 
                                            linetype = "dashed")
                        sprintf("%s threshold intensity: %s", marker, 
                                as.character(selected_valley_xcord[[marker]]))
                    }
                    else {
                        p <- p + geom_vline(aes(xintercept = log(marker_threshold)), 
                                            linetype = "dashed")
                        sprintf("%s threshold intensity: %.2f", marker, 
                                marker_threshold)
                    }
                    p <- p + theme_bw()
                    print(p)
                }
            }
        }
    }
    else {
        for (marker in markers) {
            marker_status_name <- paste(marker, "_status", sep = "")
            marker_pred_exp_colname <- paste(marker, "_predicted_phenotype", 
                                             sep = "")
            accuracy_df <- data.frame(rep(NA, nrow(predicted_data)))
            colnames(accuracy_df) <- "status"
            marker_exp_pred <- predicted_data[, c(marker_pred_exp_colname)]
            accuracy_df <- cbind(accuracy_df, marker_exp_pred)
            accuracy_df <- data.frame(accuracy_df[, "status"])
            colnames(accuracy_df) <- "status"
            marker_specific_level <- predicted_data[, marker]
            marker_specific_level <- data.frame(marker_specific_level)
            colnames(marker_specific_level) <- "Marker_level"
            level_and_accuracy <- cbind(marker_specific_level, 
                                        accuracy_df)
            if (plot_distribution) {
                p <- ggplot(level_and_accuracy, aes(x = log(Marker_level))) + 
                    geom_density()
                title <- paste("Density distribution of", marker, 
                               sep = " ")
                p <- p + labs(title = title, x = "Level of intensity (log)", 
                              y = "Density")
                if (!is.null(selected_valley_xcord[[marker]])) {
                    p <- p + geom_vline(aes(xintercept = log(selected_valley_xcord[[marker]])), 
                                        linetype = "dashed")
                    print(paste(marker, " threshold intensity: ", 
                                selected_valley_xcord[[marker]]))
                }
                else {
                    p <- p + geom_vline(aes(xintercept = log(marker_threshold)), 
                                        linetype = "dashed")
                    print(paste(marker, " threshold intensity: ", 
                                marker_threshold))
                }
                p <- p + theme_bw()
                print(p)
            }
        }
    }
    phenotype_predictions <- predicted_data[, grep("_predicted_phenotype", 
                                                   colnames(predicted_data))]
    colnames(phenotype_predictions) <- gsub("_predicted_phenotype", 
                                            "", colnames(phenotype_predictions))
    phenotype_predictions_collapsed <- vector()
    for (marker in markers) {
        temp <- phenotype_predictions[, marker]
        temp[temp == 1] <- marker
        phenotype_predictions_collapsed <- cbind(phenotype_predictions_collapsed, 
                                                 temp)
    }
    phenotype_predictions_vector <- apply(phenotype_predictions_collapsed, 
                                          1, paste, collapse = ",")
    phenotype_predictions_vector <- gsub("^0", "", phenotype_predictions_vector)
    phenotype_predictions_vector <- gsub(",0", "", phenotype_predictions_vector)
    phenotype_predictions_vector <- gsub("NA,", "", phenotype_predictions_vector)
    phenotype_predictions_vector <- gsub(",NA", "", phenotype_predictions_vector)
    phenotype_predictions_vector <- gsub("^,", "", phenotype_predictions_vector)
    phenotype_predictions_vector[phenotype_predictions_vector == 
                                     ""] <- "None"
    if (!reference_phenotypes) {
        SummarizedExperiment::colData(sce_object)$Phenotype <- phenotype_predictions_vector
        return_results <- sce_object
    }
    else {
        return_results <- predicted_data
    }
    return(return_results)
}
