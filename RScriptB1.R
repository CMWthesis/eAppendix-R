library(ape)
library(dplyr)
library(magrittr)


args <- commandArgs(TRUE)

split <- strsplit(args[1], '[.]')
Pa <- split[[1]][1]
Px <- split[[1]][2]
print(getwd())
getList <- function(file_name_pattern, ...) {
  file_list_vector  <- list.files(pattern = "bestTree") # makes a list of all files with the ending ".pestPG" which is the Delta file
  FileList <- lapply (seq_along(file_list_vector), FUN = function (x) {
    if (length(grep(file_list_vector[x], pattern = file_name_pattern)) == 1) {
      file_list_vector[x]
    }
  })
  NullFilter <- !sapply(FileList, is.null)
  FileList[NullFilter]
}
print(paste("read in", args[1], "trees", sep = " "))
tree_list <- getList("RAxML")

dMatrix_list <- tree_list %>% lapply(FUN = function(x){
  tree <- read.tree(x)
  matrix <- cophenetic(tree) %>% as.data.frame()
})


dMatrix_list2 <- tree_list %>% lapply(FUN = function(x){
  tree <- read.tree(x)})
  
  
names(dMatrix_list) <- tree_list
print(paste0("analyse", args[1], "trees"))
testLengths <- dMatrix_list %>% lapply(FUN = function(matrix){
  rd1 <- matrix[Px, Pa]
  rd2 <- matrix[Px, "18_ORa1_PxNQ"]
  rd3 <- matrix[Px, "19_MRa2_PxNQ"]

  if(rd3 > rd1 & rd2 > rd1){print('Admix')} #PxA-PxH > PxA-Pa
}) %>% do.call(rbind, .) %>% as.data.frame()
print(paste("make ratios for", args[1], "trees", sep = " "))
tRatios <- dMatrix_list %>% lapply(FUN = function(matrix){
  d2 <- matrix[Pa, "18_ORa1_PxNQ"]
  d3 <- matrix[Pa, Px]
  d4 <- matrix[Pa, "19_MRa2_PxNQ"]
  d1 <- mean(d2, d4)
  d3/d1
}) %>%
  do.call(rbind, .) %>%
  as.data.frame()

print(paste("summarized", args[1], "ratios", sep = " "))
summary <- tRatios %>% split(cut(tRatios[[1]], seq(0, 2, 0.2))) %>%
  lapply(function(x){
    nrow(x)
  }) %>%
  do.call(rbind, .) %>%
  cbind(sum(.)) %>%
  set_colnames(., c("number_in_bin", "proportion_of_total"))

for (i in 1:10) {
  summary[i,2] <- summary[i,1]/summary[i,2]
}

i <- sum(tRatios[1] >= 1)
j <- sum(tRatios[1] <= 1)
Rstar <- (i-j)/(i+j)


print(paste("Write", args[1], "files", sep = " "))
name <- paste(Pa, Px, sep = ".")
full_path_tRatios_summary <- paste0("/home/chris/Desktop/allSites_50kb/summary/", name, ".summary")
full_path_tRatios <- paste0("/home/chris/Desktop/allSites_50kb/ratios/", name, ".dist.prop")
full_path_admixTrees <- paste0("/home/chris/Desktop/allSites_50kb/admixTree_regions/", name, ".txt")
full_path_densitree <- paste0("/home/chris/Desktop/allSites_50kb/densitree/", name, ".txt")
full_path_tRatios_Rstar <- paste0("/home/chris/Desktop/allSites_50kb/summary/", name, ".Rstar")
dMatrix_list2 %>% do.call(rbind, .) %>%
  write.table(file = full_path_densitree,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

summary %>% write.table(file = full_path_tRatios_summary,
                        quote = FALSE,
                        row.names = FALSE)
Rstar %>% write.table(file = full_path_tRatios_Rstar,
                        quote = FALSE,
                        row.names = FALSE)
tRatios %>% write.table(file = full_path_tRatios,
                        quote = FALSE,
                        col.names = FALSE)

rownames(testLengths) %>% rbind() %>% t() %>% write.table(file = full_path_admixTrees,
                                                          quote = FALSE,
                                                          row.names = FALSE,
                                                          col.names = FALSE)



