library(ggplot2)
library(spatstat.random)
library(spaSim)
set.seed(610)
X <- rmpoispp(c(0.0002,0.0003,0.0005),win = owin(xrange = c(0,2000), 
                                                 yrange = c(0,2000)), 
              types = c("Tumour","Immune", "Others"))

df <- data.frame(X)
svglite::svglite("Results/mix_bg-spatstat.svg", width = 5.5, height = 4)
ggplot(df, aes(x, y, color = marks)) + geom_point(size = 1) + 
    scale_color_manual(values=c("red", "darkgreen", "lightgray"))+
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
dev.off()


length(which(df$marks == "Tumour")) # 818
length(which(df$marks == "Immune")) # 1169
length(which(df$marks == "Others")) # 1992
# total: 3979


svglite::svglite("Results/mix_bg-spaSim.svg", width = 5.5, height = 4)
set.seed(61)
Y <- TIS(n_cells = 4000, bg_method = "Hardcore",
         width = 2000, height = 2000, min_d = 10, 
         oversampling_rate = 1.5, names_of_bg_cells = c("Tumour", "Immune", "Others"), 
         proportions_of_bg_cells = c(0.2, 0.3, 0.5), plot_image = T, 
         plot_categories = c("Tumour", "Immune", "Others"),
         plot_colours = c("red", "darkgreen", "lightgray"))
dev.off()

length(which(Y$Cell.Type == "Tumour")) #779
length(which(Y$Cell.Type == "Immune")) #1209
length(which(Y$Cell.Type == "Others")) #2012
#total:4000


set.seed(1)
Z2 <- rCauchy(4, 0.01, 500)
df3 <- data.frame(Z2)
ggplot(df3, aes(x, y)) + geom_point(size = 0.5)+
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())



