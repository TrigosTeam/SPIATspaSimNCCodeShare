# this script is for simulating evenly spaced cells
library(spaSim)
load("Objects/diabetes_real_spe.Rda")
even_bg_paired <- list()

avg_d_real <- c()
for ( i in 1:length(Real)){
    image <- Real[[i]]
    avg_d_real <- c(avg_d_real, SPIAT::average_minimum_distance(image))
}

for ( i in 1:845){
    image <- Real[[i]]
    n_real <- dim(image)[2]
    x_real <- max(SpatialExperiment::spatialCoords(image)[, "Cell.X.Position"]) - 
        min(SpatialExperiment::spatialCoords(image)[, "Cell.X.Position"])
    y_real <- max(SpatialExperiment::spatialCoords(image)[, "Cell.Y.Position"]) - 
        min(SpatialExperiment::spatialCoords(image)[, "Cell.Y.Position"])
    d_real <- avg_d_real[i]
    
    LS <- 100
    for (j in seq(0.25, 0.32, 0.005)){
        image <- simulate_background_cells(n_real, x_real, y_real, method = "Even", jitter = j)
        image <- SPIAT::format_colData_to_spe(image)
        d_sim <- SPIAT::average_minimum_distance(image)
        LS_temp <- (d_sim - d_real)**2
        if (LS_temp < LS) {
            LS <- LS_temp
            j_best <- j
            print(paste("j_best", j_best))
            image_best <- image
        }
    }
    
    even_bg_paired[[i]] <- list(image = image_best, jitter  = j_best)
    # save(even_bg_paired, file = "Objects/even_bg_paired.Rda")
    
    print(i)
}



