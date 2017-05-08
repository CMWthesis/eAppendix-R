library(dplyr)
library(ggplot2)
library(gsubfn)

CHRO <- "KB2077.1"
start <- 100000
end <- 150000

PxCalca.v.PaCalca <- read.table(file = "PlutellaSA_PaCalcavPxCalca.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxCalca vs PaCalca", winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)))

  PxBB.v.PaBB <- read.table(file = "PlutellaSA_PaBBvPxBB.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxBB vs PaBB", winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)))

PxH.v.PaBB <- read.table(file = "PlutellaSA_PxHawvPaBB.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaBB", winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)))

PxH.v.PaCalca <- read.table(file = "PlutellaSA_PxHawvPaCalca.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaCalca", winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)))

PxH.v.PxBB <- read.table(file = "PlutellaSA_PxHawvPxBB.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PxBB", winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)))

  PxH.v.PxCalca <- read.table(file = "PlutellaSA_PxHawvPxCalca.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PxCalca", winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)))


cols <- c("region"="grey")
ggplot() + 
  scale_fill_manual(name="loci",values=cols) +
  geom_rect(aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = "region")) + 
  geom_line(data = PxCalca.v.PaCalca, aes(x = winCenter, y = MEAN_FST, colour = Populations), size = 1) +
  geom_line(data = PxBB.v.PaBB, aes(x = winCenter, y = MEAN_FST, colour = Populations), size = 1) +
  geom_line(data = PxH.v.PaBB, aes(x = winCenter, y = MEAN_FST, colour = Populations), size = 1) +
  geom_line(data = PxH.v.PaCalca, aes(x = winCenter, y = MEAN_FST, colour = Populations), size = 1) +
  geom_line(data = PxH.v.PxBB, aes(x = winCenter, y = MEAN_FST, colour = Populations), size = 1) +
  geom_line(data = PxH.v.PxCalca, aes(x = winCenter, y = MEAN_FST, colour = Populations), size = 1) + 
  theme_linedraw() +
  scale_colour_brewer(palette="Dark2")


