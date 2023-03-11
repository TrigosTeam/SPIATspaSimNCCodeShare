library(SPIAT)

##### all
load("../NoteN2/Objects/PCa_spe.Rda")

# R-BT
df <- data.frame(Patient = NA, image = NA, R_BC = NA)
count <- 0
for (j in c(1:16, 18:26, 28)){
    spe_list <- PCa_spe[[j]]
    for (i in 1:length(spe_list)){
        spe <- spe_list[[i]]
        r <- R_BC(spe, "AMACR", "Phenotype")
        count <- count+1
        df[count, "Patient"] <- j
        df[count, "image"] <- i
        df[count, "R_BC"] <- round(r, 2)
        print(paste(i, j))
    }
} 

save(df, file="Objects/R_BC.Rda")


##### plot #####
v1 <- c("AMACR", "CD3,CD4", "CD3,CD8")
v2 <- c("#FC8D62", "#7570B3", "#7570B3")
v3 <- c("Tumor", "Immune", "Immune")
v4 <- c("Tumor", "Immune")
v5 <- c("#FC8D62", "#7570B3")

# for (i in c(2,9,11,13,16,21,22,26,28)){
#     load(paste("~/Desktop/patients_spe/formatted_S",i,".RData", sep = ""))
# }
# clear margin
spe<-PCa_spe[[22]][[14]]
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 1.2)

spe<-PCa_spe[[22]][[6]] # Too few immune cells
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[2]][[7]] 
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[9]][[11]] # Too few immune cells
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[13]][[1]]
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)



# no clear margin
spe<-PCa_spe[[11]][[3]] # Too few tumor cells
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[22]][[12]] # Too few tumor cells
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[26]][[13]] # Too few tumor cells
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[28]][[15]] 
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[16]][[3]] 
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[21]][[8]]  # Too few cells
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

spe<-PCa_spe[[16]][[2]] 
spe$Phenotype[which(is.na(spe$Phenotype))] <- "OTHER"
spe <- define_celltypes(spe, v1, "Phenotype", v3, "Cell.Type")
plot_cell_categories(spe, v4, v5, "Cell.Type", cex = 0.5)

