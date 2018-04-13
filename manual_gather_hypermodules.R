# manual hypermodules permutations gathering

require(data.table); require(tidyr); require(dplyr); require(reshape2); 
require(compare); require(stringr)

real <- fread("Adrenocortical_carcinoma.txt")
setwd("/.mounts/labs/reimandlab/private/users/thuang/data_3/hypermodules_implementations/Adrenocortical_carcinoma/")
paths <- list.files()
random.lst <- lapply(paths,fread)