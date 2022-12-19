library(spatstat)
library(SPIAT)
library(spaSim)
load("Objects/PCa_spe.Rda")

# Define cell types
v1 <- c("AMACR", "CD3,CD4", "CD3,CD8")
v2 <- c("red", "darkgreen", "darkgreen")
v3 <- c("Tumor", "Immune", "Immune")
v4 <- c("Tumor", "Immune")
v5 <- c("red", "darkgreen")

spe <- PCa_spe[[13]][[6]]
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")

# Plot real image
svglite::svglite("Results/real_img.svg", width = 5, height = 4)
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)
dev.off()

# Plot the locations of AMACR cells
formatted_ppp <- format_spe_to_ppp(spe)
formatted_ppp_amacr <- formatted_ppp[which(formatted_ppp$marks == "AMACR")]
svglite::svglite("Results/real_img_amacr.svg", width = 5, height = 5)
plot(formatted_ppp_amacr, cex = 0.2)
dev.off()


# fit image to the model
pp <- format_spe_to_ppp(PCa_spe[[13]][[6]])
kppm(unmark(pp[which(pp$marks == "AMACR")]), clusters = "Thomas")

kappa <- as.numeric(4.817796e-07)
scale <- as.numeric(2.048089e+02)
mu <- as.numeric(583.1824)


set.seed(610)
X <- rThomas(kappa = kappa, scale = scale, mu = mu, win = owin(xrange=c(20,2666),
                                                          yrange= c(13,1993)))
svglite::svglite("thomas.svg", width = 5, height = 4)
plot(X, cex = 0.5, pch = 20)
dev.off()




kppm(unmark(pp[which(pp$marks == "AMACR")]), clusters = "Cauchy")
kappa <- as.numeric(2.128035e-07)
scale <- as.numeric(2.140172e+02)
mu <- as.numeric(1320.304)

set.seed(610)
Y <- rCauchy(kappa = kappa, scale = scale, mu = mu, win = owin(xrange=c(20,2666),
                                                               yrange= c(13,1993)))
svglite::svglite("cauchy.svg", width = 5, height = 4)
plot(Y, cex = 0.5, pch = 20)
dev.off()


# simulate an image that is similar to this one
library(spaSim)
n_cells <- 4601
# 4601 cells
average_minimum_distance(spe)
# 17.21995
set.seed(1)
bg <- simulate_background_cells(n_cells = n_cells, method = "Hardcore",
                                width = 2646, height = 1980, 
                                min_d = 17, oversampling_rate = 1.5)

svglite::svglite("ring_spaSim.svg", width = 5, height = 4)
set.seed(610)
sim <- TIS(bg_sample = bg,
           n_double_rings = 2,
           properties_of_double_rings = 
               list(I1 = list(name_of_cluster_cell = "Tumour", size = 200, shape = "Circle", 
                              centre_loc = data.frame(x = 1200, y = 800), infiltration_types =
                                  c("Immune", "Others"), infiltration_proportions = c(0.1, 0.6),
                              name_of_ring_cell = "Tumour", immune_ring_width = 250, 
                              immune_ring_infiltration_types = c("Others", "Immune"), 
                              immune_ring_infiltration_proportions = c(0.1, 0.1),
                              name_of_double_ring_cell = "Immune", double_ring_width = 120,
                              double_ring_infiltration_types = "Others", 
                              double_ring_infiltration_proportions = 0.6), 
                    I2 = list(name_of_cluster_cell = "Tumour", size = 200, shape = "Circle",
                              centre_loc = data.frame(x = 1400, y= 600), infiltration_types = 
                                  c("Immune", "Others"), infiltration_proportions = c(0.1, 0.6),
                              name_of_ring_cell = "Tumour", immune_ring_width = 250, 
                              immune_ring_infiltration_types = c("Others", "Immune"), 
                              immune_ring_infiltration_proportions = c(0.1, 0.1),
                              name_of_double_ring_cell = "Immune", double_ring_width = 100,
                              double_ring_infiltration_types = "Others", 
                              double_ring_infiltration_proportions = 0.6)),
           plot_image = T,
           plot_categories = c("Tumour", "Immune", "Others"),
           plot_colours = c("red", "darkgreen", "lightgray"))
dev.off()

svglite::svglite("ring_spaSim-amacr.svg", width = 5, height = 4)
sim <- define_celltypes(sim, c("Tumour", "Immune", "Others"), "Cell.Type",
                        c("Tumour", "Immune", "Others"), "Phenotype")
sim_ppp <- format_spe_to_ppp(sim)
sim_ppp_tumour <- sim_ppp[which(sim_ppp$marks == "Tumour")]
plot(sim_ppp_tumour, cex = 0.2)
dev.off()
