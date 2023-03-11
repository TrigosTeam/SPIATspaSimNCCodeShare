library(spaSim)

# simulate bg cells ####
svglite::svglite("Results/simulated_imgs/bg.svg", height = 4, width = 5)
set.seed(610)
bg <- simulate_background_cells(n_cells = 5000, method = "Hardcore",
                                width = 2000, height = 2000, min_d = 10,
                                oversampling_rate = 1.5)
dev.off()
# simulate mixed bg ####
svglite::svglite("Results/simulated_imgs/mix_bg.svg", height = 4, width = 5)
set.seed(610)
sim_mix <- simulate_mixing(bg_sample = bg, 
                           idents = c("Tumour", "Immune", "Others"),
                           props = c(0.2, 0.3, 0.5),
                           plot_image = TRUE,
                           plot_colours = c("#D95F02", "#7570B3", "lightgray"))
dev.off()
# simulate clusters ####
properties_of_clusters <- list(
    C1 = list(
        name_of_cluster_cell = "Tumour", 
        size = 500,
        shape = "Circle", 
        centre_loc = data.frame(x = 400, y = 1000), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05)
    ),
    C2 = list(
        name_of_cluster_cell = "Tumour", 
        size = 650,
        shape = "Oval", 
        centre_loc = data.frame(x = 500, y = 200), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05)
    ),
    C3 = list(
        name_of_cluster_cell = "Tumour", 
        size = 500,
        shape = "Oval", 
        centre_loc = data.frame(x = 1500, y = 1600), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05)
    ),
    C4 = list(
        name_of_cluster_cell = "Immune", 
        size = 600, shape = "Irregular", 
        centre_loc = data.frame(x = 1200, y = 500), 
        infiltration_types = "Others", infiltration_proportions = 0.2)
)

svglite::svglite("Results/simulated_imgs/clusters.svg", height = 4, width = 5)
set.seed(610)
sim_cluster <- simulate_clusters(
    bg_sample = bg1,
    n_clusters = 4,
    cluster_properties = properties_of_clusters,
    plot_image = T,
    plot_categories = c("Tumour", "Immune", "Others"),
    plot_colours = c("#D95F02", "#7570B3", "lightgray"))
dev.off()


# simulate immune rings ####
properties_of_rings <- list(
    I1 = list(
        name_of_cluster_cell = "Tumour", 
        size = 550, 
        shape = "Oval", 
        centre_loc = data.frame(x = 900, y = 800), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.15, 0.1),
        name_of_ring_cell = "Immune", 
        immune_ring_width = 120,
        immune_ring_infiltration_types = "Others", 
        immune_ring_infiltration_proportions = 0.2
    ),
    I2 = list(
        name_of_cluster_cell = "Tumour", 
        size = 650,
        shape = "Circle", 
        centre_loc = data.frame(x = 850, y = 400), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.15, 0.1),
        name_of_ring_cell = "Immune", 
        immune_ring_width = 120,
        immune_ring_infiltration_types = "Others", 
        immune_ring_infiltration_proportions = 0.2
    ),
    I3 = list(
        name_of_cluster_cell = "Tumour", 
        size = 200,
        shape = "Oval", 
        centre_loc = data.frame(x = 1800, y = 900), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.15, 0.1),
        name_of_ring_cell = "Immune", 
        immune_ring_width = 85,
        immune_ring_infiltration_types = "Others", 
        immune_ring_infiltration_proportions = 0.15
    )
)

svglite::svglite("Results/simulated_imgs/rings.svg", height = 4, width = 5)
set.seed(610)
sim_ring <- simulate_immune_rings(
    bg_sample = bg1,
    n_ir = 3,
    ir_properties = properties_of_rings,
    plot_image = T,
    plot_categories = c("Tumour", "Immune", "Others"),
    plot_colours = c("#D95F02", "#7570B3", "lightgray"))
dev.off()

# simulate vessles ####
# simulate big tumor clusters first
properties_of_bigC <- list(
    C1 = list(
        name_of_cluster_cell = "Tumour", 
        size = 800,
        shape = "Circle", 
        centre_loc = data.frame(x = 400, y = 1200), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05)
    ),
    C2 = list(
        name_of_cluster_cell = "Tumour", 
        size = 850,
        shape = "Oval", 
        centre_loc = data.frame(x = 400, y = 250), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05)
    ),
    C3 = list(
        name_of_cluster_cell = "Tumour", 
        size = 500,
        shape = "Oval", 
        centre_loc = data.frame(x = 1600, y = 1700), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05)
    ),
    C4 = list(
        name_of_cluster_cell = "Tumour", 
        size = 600, shape = "Circle", 
        centre_loc = data.frame(x = 1750, y = 700), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05)),
    C5 = list(
        name_of_cluster_cell = "Tumour", 
        size = 400, shape = "Circle", 
        centre_loc = data.frame(x = 1750, y = 200), 
        infiltration_types = c("Immune", "Others"), 
        infiltration_proportions = c(0.1, 0.05))
)

set.seed(610)
sim_bigC <- simulate_clusters(
    bg_sample = bg1,
    n_clusters = 5,
    cluster_properties = properties_of_bigC,
    plot_image = T,
    plot_categories = c("Tumour", "Immune", "Others"),
    plot_colours = c("#D95F02", "#7570B3", "lightgray"))

# simulate vessels
stripe_properties = list(
    S1 = list(
        number_of_stripes = 2, 
        name_of_stripe_cell = "Others", 
        width_of_stripe = 60, 
        infiltration_types = c("Immune"),
        infiltration_proportions = c(0.08)
        )
    )

svglite::svglite("Results/simulated_imgs/vessels.svg", height = 4, width = 5)
set.seed(10)
sim_vessle <- simulate_stripes(bg_sample = sim_bigC,
                               n_stripe_type = 1,
                               stripe_properties = stripe_properties,
                               plot_image = T,
                               plot_categories = c("Tumour", "Immune", "Others"),
                               plot_colours = c("#D95F02", "#7570B3", "lightgray"))
dev.off()
