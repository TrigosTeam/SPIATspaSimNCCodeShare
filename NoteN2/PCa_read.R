# This script is for generating PCa spe objects
library(SPIAT)

dir <- "~/Documents/data/" # change this for your own path (absolute)
setwd(dir)
markers <- c("DAPI", "CD3", "PDL-1", "FOXP3", "CD4", "CD8", "AMACR")
locations <- c("Nucleus", "Cytoplasm", "Membrane", "Cytoplasm", 
               "Cytoplasm", "Cytoplasm", "Cytoplasm")

# save the objects as a list of lists
PCa_spe <- list()
for (i in c(1:16, 18:26, 28)){ # Patient level
    PCa_spe[[i]] <- list() 
    files <- list.files(paste0(dir, "SPORADIC/S", i, "/"))
    files <- files[-grep("summary",files)]
    
    for (j in seq_len(length(files))){ # Image level
        path <- paste0(dir, "SPORADIC/S", i, "/", files[j])
        spe <- format_inform_to_spe(path = path, markers = markers,
                                    locations = locations)
        attr(spe, "name") <- paste0("S", i, ".", j)
        PCa_spe[[i]][[j]] <- spe
        save(PCa_spe, file = paste0(dir, "PCa_spe.Rda"))
    }
}



