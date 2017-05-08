library(dplyr)
library(tidyr)
library(plyr)



#read in vcf
vcf <- read.table("Calca6_wolbachia_forR.sites.vcf", header = TRUE) %>% separate(col = Calca6, into = c("GT", "DP"), sep = ":")
#get the list of scaffold in the vcf
contig_list_in_vcf <- as.data.frame(vcf[1]) %>% distinct()


# separate sequence data into contigs
separate_contigs <- as.matrix(contig_list_in_vcf) %>% lapply(function(x){
  filter(vcf, CHROM == x )
})


## now we need to split the contigs into markers
## l is the length of the dataframe -1 as this is the mt scaffold

split_points <- separate_contigs %>% lapply(function(df = .){
  l <- 1:(nrow(df)-1)
  lapply(l, function(x = l, y = df){
    y[x,2] >= y[(x+1),2]-100}) %>% do.call(rbind, .) %>% rbind(1)})

make_float_points <- split_points %>% lapply(function(x){
  length <- seq_along(x)
  float <- 1
  data <- data.frame(1)
  for (i in length){
    if ((x[i] == 1)) {
      data[i,] <- float
    }
    else {
      float <- float+1
      data[i,] <- float}
  }
  data
})

split_data <- seq_along(make_float_points) %>% lapply(function(x){
  cbind(separate_contigs[[x]], make_float_points[[x]]) %>% 
    group_by(X1) %>% 
    split(.$X1)})

fasta_sequence <- split_data %>% lapply(function(x){
  x %>% lapply(function(y){
    length <- 1:nrow(y)
    fasta <- c()
    for (i in length) {
      if (!(y[i, 12] == 1)) {
        fasta[i] <- as.matrix(y[i, 4])
      }
      else {fasta[i] <- as.matrix(y[i, 5])}
    }
   data.frame(fasta) %>% t() %>% as.data.frame()
  })  %>%  do.call(rbind.fill, .)
}) %>%  do.call(rbind.fill, .)





