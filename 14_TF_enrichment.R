# transcription factor enrichment in module genes 

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(stringr)

# files
tcga.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/modules_added_nonlinkers_permut_10.rds")
network <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv", header = F)
tfs <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/Human_TF_MasterList_v1_02.csv")


# file names

FormatTFs <- function(tfs) {
	tfs$is_tf <- tfs$'Is TF?'
	filtered <- filter(tfs, is_tf == "Yes")
	tf.genes <- toupper(unique(filtered$V2))

	return(tf.genes)
}

FormatModuleGenes <- function(tcga.modules) {
	mods <- tcga.modules %>% 
				ungroup() %>%
				mutate(Module = strsplit(Module, split = ";")) %>% 
				unnest()

	mod.genes <- unique(mods$Module)

	return(mod.genes)
}

FormatBiogridGenes <- function(network) {
	net.genes <- unique(unlist(network))

	return(net.genes)
}

TFEnrichment <- function(tf.genes, mod.genes, biogrid.genes) {
	fisher <- fisher.test(biogrid.genes %in% tf.genes, biogrid.genes %in% mod.genes, alternative="greater")
	table(biogrid.genes %in% tf.genes, biogrid.genes %in% mod.genes)

	return(fisher)
}

GetExpectedTFEnrichment <- function(tf.genes, mod.genes, biogrid.genes, permut=10000) {
	sample.size <- length(mod.genes)
	sample.tf.number <- length(intersect(mod.genes, tf.genes))
	background.tf.numbers <- replicate(permut, sum(tf.genes %in% sample(biogrid.genes, sample.size, replace = TRUE)))
	background.greater <- table(background.tf.numbers >= sample.tf.number)
	background.average <- mean(background.tf.numbers)

	return(list(background.average, sample.tf.number, background.greater))
}

############################################################################################################
tf.genes <- FormatTFs(tfs)
mod.genes <- FormatModuleGenes(tcga.modules)
biogrid.genes <- FormatBiogridGenes(network)
fisher <- TFEnrichment(tf.genes, mod.genes, biogrid.genes)
tf.permutations <- GetExpectedTFEnrichment(tf.genes, mod.genes, biogrid.genes, permut=10000) # no enrichment