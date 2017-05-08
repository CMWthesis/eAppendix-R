library(ggplot2)
library(dplyr)
library(magrittr)
library(grid)
library(gtable)
library(reshape)
library(RColorBrewer)

setwd("~/Desktop/allSites_50kb/ratios")
getList <- function(file_name_pattern, ...) {
  file_list_vector  <- list.files(pattern = "prop") # makes a list of all files with the ending ".pestPG" which is the Delta file
  FileList <- lapply (seq_along(file_list_vector), FUN = function (x) {
    if (length(grep(file_list_vector[x], pattern = file_name_pattern)) == 1) {
      file_list_vector[x]
    }
  })
  NullFilter <- !sapply(FileList, is.null)
  FileList[NullFilter]
  
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##################
#### GET FILES
######################
files <- getList(file_name_pattern = "dist") %>% lapply(function(x){
  if (grepl("Cook", x)) {
    region <- "ACT_Cook_2014"
  } 
  if (grepl("GF", x)) {
    region <- "ACT_GF_2015"
  } 
  if (grepl("Calca", x)) {
    region <- "SA_2014"
  } 
  if (grepl("BB", x)) {
    region <- "SA_2014"
  } 
  read.table(x) %>% cbind(name = x, region) %>% 
    set_names(c("site", "Ratio", "Comparison", "Region"))
}) %>% do.call(rbind, .)

#########################
########## histogram
#######################
p_histogram <- ggplot(files) +
  geom_histogram(aes(x = Ratio, group = Comparison), binwidth = 0.2, alpha = 0.5) +
  xlim(0, 2) +
  geom_vline(aes(xintercept = 1))+
  geom_hline(aes(yintercept = 0)) +
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("") + 
  facet_wrap(~ Region, nrow = 1)
############################
########### density plot
#########################
p_density1 <- ggplot(files) +
  geom_density(aes(x = Ratio, colour = Comparison), size = 0.4) +
  xlim(0, 2) +
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 1)) +
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  ylab("") + 
  facet_wrap( ~ Region, nrow = 1)

setwd("~/Dropbox/ChrisWard_Honours/FINAL_pipeline_scripts/Fst")
CHRO <- "KB207380.1"
start <- 6000
end <- 7000
###############################
############# line graph
############################
scaffold <- "KB207380.1"
line_locus <- files %>% 
  filter(grepl(scaffold, .$site))
line_locus$site <- gsub("^.*_","", line_locus$site)
line_locus$site <- gsub(".phy","", line_locus$site)

KB207380line <- line_locus %>%
  ggplot() +
  scale_fill_manual(name="loci",values=cols) +
  geom_line(aes(y = Ratio, x = as.numeric(site)/100, group = Comparison)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none") +
  theme_linedraw() + 
  geom_hline(yintercept = 1) +
  ylim(0, 2) +
  geom_rect(aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = "100kb region")) + 
  facet_wrap( ~ Region, nrow = 1) 

########
#fst plot
########
setwd("~/Dropbox/ChrisWard_Honours/FINAL_pipeline_scripts/Fst")
CHRO <- "KB207758.1"
start <- 6000
end <- 7000

PxSA.v.PxH <- read.table(file = "PlutellaSA_PxSA.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaH", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "SA 2014")

PxSA.v.PaSA <- read.table(file = "PlutellaSA_PxSA.PaSA.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "SA 2014")

PxH.v.PaSA <- read.table(file = "PlutellaSA_PaSA.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "SA 2014")

PxACT14.v.PxH <- read.table(file = "PlutellaACT14_PxCook.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaH",
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2014")

PxACT14.v.PaACT14 <- read.table(file = "PlutellaACT14_PxCook.PaCook.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaA",
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2014")

PxH.v.PaACT14 <- read.table(file = "PlutellaACT14_PaCook.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2014")

PxACT15.v.PxH <- read.table(file = "PlutellaACT15_PxGF.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaH",
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2015")

PxACT15.v.PaACT15 <- read.table(file = "PlutellaACT15_PxGF.PaGF.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2015")

PxH.v.PaACT15 <- read.table(file = "PlutellaACT15_PaSA.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaA", winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE))) %>%
  cbind(demulti = "ACT 2015")

merged <- rbind(PxSA.v.PxH, PxSA.v.PaSA, PxH.v.PaSA,
                PxACT14.v.PxH, PxACT14.v.PaACT14, PxH.v.PaACT14,
                PxACT15.v.PxH, PxACT15.v.PaACT15, PxH.v.PaACT15)

cols <- c("100kb region"="grey")

KB207380fst <- ggplot(data = merged) + 
  scale_fill_manual(name="loci",values=cols) +
  geom_rect(aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = "100kb region")) + 
  geom_line( aes(x = as.numeric(winCenter)/100, y = WEIGHTED_FST, 
                               colour = Populations), size = 1) +
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none") +
  facet_wrap( ~ demulti, nrow = 1)


############################
scaffold <- "KB207380.1"
line_locus <- files %>% 
  filter(grepl(scaffold, .$site))
line_locus$site <- gsub("^.*_","", line_locus$site)
line_locus$site <- gsub(".phy","", line_locus$site)

KB207380line <- line_locus %>%
  ggplot() +
  scale_fill_manual(name="loci",values=cols) +
  geom_rect(aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = "100kb region")) + 
  theme_linedraw() + 
  geom_line(aes(y = Ratio, x = as.numeric(site)/100, group = Comparison)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none") +
  geom_hline(yintercept = 1) +
  ylim(0, 2) +
  facet_wrap( ~ Region, nrow = 1) 

########
#fst plot
########
setwd("~/Dropbox/ChrisWard_Honours/FINAL_pipeline_scripts/Fst")
CHRO <- "KB207758.1"
start <- 1441.65
end <- 1621.69

PxSA.v.PxH1 <- read.table(file = "PlutellaSA_PxSA.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaH", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "SA 2014")

PxSA.v.PaSA1 <- read.table(file = "PlutellaSA_PxSA.PaSA.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "SA 2014")

PxH.v.PaSA1 <- read.table(file = "PlutellaSA_PaSA.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "SA 2014")

PxACT14.v.PxH1 <- read.table(file = "PlutellaACT14_PxCook.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaH",
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2014")

PxACT14.v.PaACT141 <- read.table(file = "PlutellaACT14_PxCook.PaCook.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaA",
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2014")

PxH.v.PaACT141 <- read.table(file = "PlutellaACT14_PaCook.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2014")

PxACT15.v.PxH1 <- read.table(file = "PlutellaACT15_PxGF.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaH",
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2015")

PxACT15.v.PaACT151 <- read.table(file = "PlutellaACT15_PxGF.PaGF.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>%
  cbind(., Populations = "PxA vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2015")

PxH.v.PaACT151 <- read.table(file = "PlutellaACT15_PaSA.PxHaw.windowed.weir.fst", header = TRUE) %>% 
  filter(CHROM == CHRO) %>% 
  cbind(., Populations = "PxH vs PaA", 
        winCenter = rowMeans(subset(., select = c(BIN_START, BIN_END), na.rm = TRUE)), demulti = "ACT 2015")

merged <- rbind(PxSA.v.PxH1, PxSA.v.PaSA1, PxH.v.PaSA1,
                PxACT14.v.PxH1, PxACT14.v.PaACT141, PxH.v.PaACT141,
                PxACT15.v.PxH1, PxACT15.v.PaACT151, PxH.v.PaACT151)

cols <- c("100kb region"="grey")

KB207758fst <- ggplot(data = merged) + 
  scale_fill_manual(name="loci",values=cols) +
  geom_rect(aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = "100kb region")) + 
  geom_line( aes(x = as.numeric(winCenter)/100, y = WEIGHTED_FST, 
                 colour = Populations), size = 1) +
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none") +
  facet_wrap( ~ demulti, nrow = 1)

scaffold1 <- "KB207758.1"
line_locus <- files %>% 
  filter(grepl(scaffold1, .$site))
line_locus$site <- gsub("^.*_","", line_locus$site)
line_locus$site <- gsub(".phy","", line_locus$site)

KB207758line <- line_locus %>% ggplot() +
  scale_fill_manual(name="loci",values=cols) +
  geom_line(aes(y = Ratio, x = as.numeric(site)/100, group = Comparison)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none") +
  geom_hline(yintercept = 1) +
  ylim(0, 2) +
  geom_rect(aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill = "100kb region")) + 
  facet_wrap( ~ Region, nrow = 1) 



KB207380fst
KB207380line

multiplot(KB207380line, KB207380fst)
multiplot(p_density1, p_histogram)
