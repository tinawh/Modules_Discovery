# pathway enrichment 

require(data.table); require(tidyr); require(dplyr); require(reshape2)
require(gProfileR)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-19-18/"
dir.create(loc)

# files
tcga.modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/modules_added_nonlinkers_permut_10.rds")
network <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv", header = F)

# file names
total.pathway.enrichment.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-19-18/gprofiler_all_modules_PE.txt"

GetModuleGenes <- function(tcga.modules) {
	# get vector of total module genes

	unnested <- tcga.modules %>% 
					ungroup() %>%
					mutate(Module = strsplit(Module, split = ";")) %>% 
					unnest() 

	mod.genes <- unique(unlist(unnested$Module))

	return(mod.genes)
}

GetNetworkGenes <- function(network) {
	# get vector of all genes in biogrid network 

	network.genes <- unique(unlist(network))

	return(network.genes)
}

WritePathways <- function(mod.genes, network.genes, file.name = total.pathway.enrichment.file) {
	# obtain gprofiler pathway enrichment and write to file.name

	pathways <- gprofiler(mod.genes, custom_bg = network.genes, src_filter = c("GO:BP", "REAC"), ordered_query = TRUE, exclude_iea = TRUE, min_isect_size = 10, min_set_size = 10, max_set_size = 300)
	subsets <- pathways %>% 
					# mutate(amount.overlap = overlap.size/term.size) %>% 
					# filter(amount.overlap > .4) %>%
					mutate(FDR = p.value) %>%
					select(term.id, term.name, p.value, FDR, query.number, intersection) %>% 
					rename(GO.ID = term.id, Description = term.name, p.Val = p.value, Phenotype = query.number, Genes = intersection)
	
	write.table(subsets, file = total.pathway.enrichment.file, quote = F, row.names = F, sep = "\t")
}
########################################################################################################################
mod.genes <- GetModuleGenes(tcga.modules)
network.genes <- GetNetworkGenes(network)
WritePathways(mod.genes, network.genes)