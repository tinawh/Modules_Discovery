# kinases enrichment in module genes 

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(stringr)

# files
tcga.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/modules_added_nonlinkers_permut_10.rds")
network <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv", header = F)
kinome <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/kinome_known_kinases.txt")

# file names

FormatKinomeGenes <- function(kinome) {
	kinome.kinases <- toupper(unique(kinome$Name))

	return(kinome.kinases)
}

FormatModuleGenes <- function(tcga.modules) {
	mods <- tcga.modules %>% 
				ungroup() %>%
				mutate(Module = strsplit(Module, split = ";")) %>% 
				unnest()

	mod.genes <- toupper(unique(mods$Module))

	return(mod.genes)
}

FormatBiogridGenes <- function(network) {
	net.genes <- unique(unlist(network))

	return(net.genes)
}

KinaseEnrichment <- function(kinome.kinases, mod.genes, biogrid.genes) {
	fisher <- fisher.test(biogrid.genes %in% kinome.kinases, biogrid.genes %in% mod.genes, alternative="greater")
	table(biogrid.genes %in% kinome.kinases, biogrid.genes %in% mod.genes)

	return(fisher)
}

GetExpectedKinaseEnrichment <- function(kinome.kinases, mod.genes, biogrid.genes, permut=10000) {
	sample.size <- length(mod.genes)
	sample.kinases.number <- length(intersect(mod.genes, kinome.kinases))
	background.kinases.numbers <- replicate(permut, sum(kinome.kinases %in% sample(biogrid.genes, sample.size, replace = TRUE)))
	background.greater <- table(background.kinases.numbers >= sample.kinases.number)
	background.average <- mean(background.kinases.numbers)

	return(list(background.average, sample.kinases.number, background.greater))
}

############################################################################################################
kinome.kinases <- FormatKinomeGenes(kinome)
mod.genes <- FormatModuleGenes(tcga.modules)
biogrid.genes <- FormatBiogridGenes(network)
fisher <- KinaseEnrichment(kinome.kinases, mod.genes, biogrid.genes)
kinome.permutations <- GetExpectedKinaseEnrichment(kinome.kinases, mod.genes, biogrid.genes, permut=10000)
