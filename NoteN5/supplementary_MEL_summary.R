load("../Figure6/Objects/Pan_immune_spe.RData")
library(SPIAT)

plot_cell_categories(primary, c("Undefined","MEL", "T", "B", "DC"),
                     c("lightgray", "red", "blue", "darkgreen", "orange"),
                     feature_colname = "Cell.Type", cex = 0.4, layered = TRUE)

plot_cell_categories(MET, c("Undefined","MEL", "T", "B", "DC"),
                     c("lightgray", "red", "blue", "darkgreen", "orange"),
                     feature_colname = "Cell.Type", cex = 0.2, layered = TRUE)

plot_cell_categories(REL, c("Undefined","MEL", "T", "B", "DC"),
                     c("lightgray", "red", "blue", "darkgreen", "orange"),
                     feature_colname = "Cell.Type", cex = 0.2, layered = TRUE)

# Code for cell proportions see "../Figure6/Fig6-Melanoma_margin_detection.Rmd"