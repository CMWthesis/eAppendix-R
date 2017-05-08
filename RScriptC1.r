library(ape)
library(dplyr)
library(magrittr)
library(ggplot2)
library(gplots)

culex <- read.csv("culex_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "CDS" & !grepl("hypothetical", V9))
culex$V9 <- gsub("^.*=","", culex$V9)
culex_names <- unique(culex[[9]])

australiana <- read.csv("australiana_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "CDS" & !grepl("hypothetical", V9))
australiana$V9 <- gsub("^.*=","", australiana$V9)
australiana_names <- unique(australiana[[9]])

sim <- read.csv("simulans_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "CDS" & !grepl("hypothetical", V9))
sim$V9 <- gsub("^.*=","", sim$V9)
sim_names <- unique(sim[[9]])

melano <- read.csv("melanogaster_annotation.csv", fill = TRUE, comment.char = "", header = FALSE) %>%
  filter(V3 == "CDS" & !grepl("hypothetical", V9))
melano$V9 <- gsub("^.*=","", melano$V9)
melano_names <- unique(melano[[9]])

merged_names <- c(culex_names, melano_names, sim_names) %>% unique() 

aus_in <- merged_names %in% australiana_names %>% as.numeric() %>% set_names(merged_names)
culex_in <- merged_names %in% culex_names %>% as.numeric() %>% set_names(merged_names)
sim_in <- merged_names %in% sim_names %>% as.numeric() %>% set_names(merged_names)
melano_in <- merged_names %in% melano_names %>% as.numeric() %>% set_names(merged_names) 


sim_df <- data.frame(number = sim_in, delim = "wSim", gene = seq_along(merged_names))
mel_df <- data.frame(number = melano_in, delim = "wMel", gene = seq_along(merged_names))
aus_df <- data.frame(number = aus_in, delim = "wAus", gene = seq_along(merged_names))
culex_df <- data.frame(number = culex_in, delim = "wPip", gene = seq_along(merged_names))

names_df <- rbind(sim_df, mel_df, aus_df, culex_df)

names_df %>%
  ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(aes(x = gene, y = number, colour = delim), alpha = 0.5) +
  facet_wrap( ~ delim, nrow = 4)

names_df %>%
  unique() %>%
  ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_bar(aes(x = gene, y = number), fill = "black", colour = "black", stat = "identity", size = 0.3) +
  facet_wrap( ~ delim, nrow = 4)



my_palette <- colorRampPalette(c("white", "black", "black"))(n = 1)
heatmap.2(names_df, Colv = NA, black(10), 
          density.info="none", trace="none", 
          symm=F,symkey=F,symbreaks=F, scale="none", na.color = "black")


