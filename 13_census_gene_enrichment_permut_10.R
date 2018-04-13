# cancer census enrichment in module genes 

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(stringr)

# files
tcga.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/modules_added_nonlinkers_permut_10.rds")
network <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv", header = F)
census <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/Census_allMon Jan 22 18_25_23 2018.tsv")

# file names

FormatCensusGenes <- function(census) {
	census.genes <- toupper(unique(census$"Gene Symbol"))

	return(census.genes)
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

CensusEnrichment <- function(census.genes, mod.genes, biogrid.genes) {
	fisher <- fisher.test(biogrid.genes %in% census.genes, biogrid.genes %in% mod.genes, alternative="greater")
	table(biogrid.genes %in% census.genes, biogrid.genes %in% mod.genes)

	return(fisher)
}

GetExpectedCensusEnrichment <- function(census.genes, mod.genes, biogrid.genes, permut=10000) {
	sample.size <- length(mod.genes)
	sample.census.number <- length(intersect(mod.genes, census.genes))
	background.census.numbers <- replicate(permut, sum(census.genes %in% sample(biogrid.genes, sample.size, replace = TRUE)))
	background.greater <- table(background.census.numbers >= sample.census.number)
	background.average <- mean(background.census.numbers)

	return(list(background.average, sample.census.number, background.greater))
}

############################################################################################################
census.genes <- FormatCensusGenes(census)
mod.genes <- FormatModuleGenes(tcga.modules)
biogrid.genes <- FormatBiogridGenes(network)
fisher <- CensusEnrichment(census.genes, mod.genes, biogrid.genes)

census.permutations <- GetExpectedCensusEnrichment(census.genes, mod.genes, biogrid.genes, permut=10000)
