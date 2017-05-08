library(rgl)
library(dplyr)
library(magrittr)
library(RColorBrewer)
library(pca3d)
library(ggplot2)


EigenVecs <- read.table("1-29.WGS.evec") %>% extract(1:4) %>% set_colnames(c("Collection Site", "EV1", "EV2", "EV3"))

plotData1 <- select(EigenVecs, Ind, EV1 = EV1, EV2 = EV2) %>%
  cbind(., PCA_num = "EigV1 vs EigV2")
plotData2 <- select(EigenVecs, Ind, EV1 = EV1, EV2 = EV3)  %>%
  cbind(., PCA_num = "EigV1 vs EigV3")
plotData3 <- select(EigenVecs, Ind, EV1 = EV2, EV2 = EV3) %>%
  cbind(., PCA_num = "EigV2 vs EigV3")

plotData_Final <- rbind(plotData1)

plotData_Final %>% ggplot(aes(x = EV1, y = EV2, colour = Ind )) +
  geom_point(size = 3) +
  theme_bw() +
  facet_wrap(~PCA_num) +
  xlab("Eigenvalue x") +
  ylab("Eigenvalue y")


EigenVecs[2:4] %>%
  as.matrix() %>%
  pca3d(group = as.matrix(EigenVecs[1]),
        axes.color= "black",
        show.scale = TRUE,
        radius = 2,
        show.labels = FALSE, new = TRUE,
        fancy = FALSE, 
        show.shapes = FALSE)

dev.off()



