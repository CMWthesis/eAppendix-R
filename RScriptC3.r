library(ape)
library(dplyr)
library(magrittr)
library(ggplot2)
library(gplots)

culex <- read.csv("culex_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "tRNA" & !grepl("hypothetical", V9))
culex$V9 <- gsub("^.*=","", culex$V9)
culex_names <- unique(culex[[9]])


australiana <- read.csv("australiana_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "tRNA" & !grepl("hypothetical", V9))
australiana$V9 <- gsub("^.*=","", australiana$V9)
australiana_names <- unique(australiana[[9]])

sim <- read.csv("simulans_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "tRNA" & !grepl("hypothetical", V9))
sim$V9 <- gsub("^.*=","", sim$V9)
sim_names <- unique(sim[[9]])

melano <- read.csv("melanogaster_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "tRNA" & !grepl("hypothetical", V9))
melano$V9 <- gsub("^.*=","", melano$V9)
melano_names <- unique(melano[[9]])

merged_names <- c(culex_names, melano_names, sim_names) %>% unique() 

aus_in <- merged_names %in% australiana_names %>% as.numeric() %>% set_names(merged_names)
culex_in <- merged_names %in% culex_names %>% as.numeric() %>% set_names(merged_names)
sim_in <- merged_names %in% sim_names %>% as.numeric() %>% set_names(merged_names)
melano_in <- merged_names %in% melano_names %>% as.numeric() %>% set_names(merged_names) 


names_df <- data.frame(sim_in, aus_in, culex_in,  melano_in, malayi_in) %>% t() %>% as.matrix()

my_palette <- colorRampPalette(c("green", "black", "blue"))(n = 3)
heatmap.2(names_df, Colv = NA, col=redblue(75), 
          density.info="none", trace="none", 
          symm=F,symkey=F,symbreaks=F, scale="none", na.color = "black")


