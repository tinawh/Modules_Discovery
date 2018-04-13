require(data.table); require(tidyr); require(dplyr); require(reshape2); require(compare)
require(stringr); require(igraph)

# files
network <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv", header = F)
modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/modules_added_nonlinkers_permut_10.rds")

kinome <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/kinome_known_kinases.txt")
census <- fread("/.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/Census_allMon Jan 22 18_25_23 2018.tsv")
muts <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/unsplit_mutation_data.rds")

# file names 
formatted.muts.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/formatted_muts_for_igraph.rds"
links.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/links.rds"
nodes.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/nodes.rds"
module.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/module_networks.pdf"

FormatMutationsbyCancerType <- function(muts, modules) {
	splitted.mods <- split(modules, f=list(modules$cancer, modules$Module), sep="~", drop=TRUE)
	
	mod.mutations <- list()
	for (mod in names(splitted.mods)) {
		mod.genes <- unlist(strsplit(unique(splitted.mods[[mod]]$Module), split = ";"))

		# unnest patients
		unnested.patients <- splitted.mods[[mod]] %>% 
			mutate(patient_count = strsplit(List_of_patients, split = ",")) %>% 
			unnest(patient_count) 

 		joined <- inner_join(unnested.patients, muts, by = c("cancer" = "Study_Name", "patient_count" = "patient"))

 		filtered <- joined %>% 
 			filter(hugo_name %in% mod.genes) %>% 
 			group_by(hugo_name) %>% 
 			summarize(patient_count = n()) %>% 
 			full_join(data.frame(mod.genes, stringsAsFactors=F), by = c("hugo_name" = "mod.genes")) %>% 
			mutate(patient_count = ifelse(is.na(patient_count), 1, patient_count + 1))

		mod.mutations[[mod]] <- filtered
	}
	return(mod.mutations)
}

FormatKinomeGenes <- function(kinome) {
	kinome.kinases <- toupper(unique(kinome$Name))
	
	return(kinome.kinases)
}

FormatCensusGenes <- function(census) {
	census.genes <- toupper(unique(census$"Gene Symbol"))

	return(census.genes)
}

GetiGraphNetworkLinks <- function(nonlinkers.added, network) {
	# get links 

	splitted <- split(nonlinkers.added, f = nonlinkers.added$cancer)
	arranged <- lapply(splitted, function(dframe) dframe %>% arrange(coxph))

	sorted <- bind_rows(arranged) # order 

	selected <- sorted %>%
					select(Module, cancer, intersect) %>% 
					mutate(mod.genes = strsplit(Module, split = ";")) %>% 
					unnest()

	joined <- inner_join(selected, network, by = c("mod.genes" = "V1"))

	# split by module 
	splitted.2 <- split(joined, f = joined$Module, drop = TRUE)

	# reorder
	ex <- splitted.2[sorted$Module]

	links <- list()
	for (mod in names(ex)) {

		cancer <- unique(ex[[mod]]$cancer)
		module <- unique(ex[[mod]]$Module)
		id <- paste(cancer, module, sep = "~")

		# links
		mod.dframe <- data.frame(genes = unlist(strsplit(mod, split = ";")), stringsAsFactors = FALSE)

		# remove non module genes in substrates 
		joined <- inner_join(ex[[mod]][, c("mod.genes","V2")], mod.dframe, by = c("V2" = "genes")) %>% 
						filter(mod.genes != V2)

		# remove duplication
		new.dframe <- as.data.frame(unique(t(apply(joined, 1, sort))), stringsAsFactors = FALSE)

		# rename 
		renamed.new.dframe <- rename(new.dframe, from = mod.genes, to = V2)

		links[[id]] <- renamed.new.dframe

	}

	return(links)
}

GetiGraphNetworkNodes <- function(nonlinkers.added = modules, network, links, census.genes, kinome.genes, formatted.muts) {
	# get links 

	splitted <- split(nonlinkers.added, f = nonlinkers.added$cancer)
	arranged <- lapply(splitted, function(dframe) dframe %>% arrange(coxph))

	sorted <- bind_rows(arranged) # order 

	selected <- sorted %>%
					select(Module, cancer, intersect) %>% 
					mutate(mod.genes = strsplit(Module, split = ";")) %>% 
					unnest()

	joined <- inner_join(selected, network, by = c("mod.genes" = "V1"))

	# split by module 
	splitted.2 <- split(joined, f = list(joined$cancer, joined$Module), sep = "~", drop = TRUE)

	# reorder
	ex <- splitted.2[paste(sorted$cancer, sorted$Module, sep = "~")]

	nodes <- list()
	for (dframe in names(links)) {
		intersect <- unlist(strsplit(unique(ex[[dframe]]$intersect), split = ";"))
		unlisted <- data.frame(mod.genes = unique(unlist(links[[dframe]])), stringsAsFactors=F)
		cancer <- gsub("~.*", "", dframe)

		# added is_linker, is_kinase, is_cancer_gene column
		added.gene.info <- mutate(unlisted, 
			is_linker = as.integer(!mod.genes %in% intersect) + 1, 
			is_kinase = as.integer(mod.genes %in% kinome.genes) + 1, 
			is_cancer_gene = as.integer(mod.genes %in% census.genes) + 1)

		# add number of mutations
		added.mutations <- inner_join(added.gene.info, formatted.muts[[dframe]], by = c("mod.genes" = "hugo_name"))

		# scale mutation numbers 
		scaled.mutations <- mutate(added.mutations, patient_count = ((log2(patient_count)+1)*4))

		nodes[[dframe]] <- scaled.mutations
	}
	return(nodes)
}

DrawiGraphModules <- function(links, nodes, modules=module.file) {
	pdf(module.file, width = 35, height = 36, onefile=T)

	for (mod in names(links)) {
		net <- graph_from_data_frame(d=links[[mod]], vertices=nodes[[mod]], directed=F)
		colrs <- c("#ff6666", "gray80")
		shapes <- c("circle", "square")
		label.colrs <- c("gray40", "#4f7cb3")
		V(net)$color <- colrs[V(net)$is_linker]
		V(net)$size <- V(net)$patient_count+1
		V(net)$shape <- shapes[V(net)$is_kinase]
		V(net)$label.color <- label.colrs[V(net)$is_cancer_gene]
		V(net)$label.cex <- 3
		V(net)$vertex.frame.size <- 2
		V(net)$label.family <- "sans"
		E(net)$width <- 2

		par(family = "sans")

		l <- layout_with_kk(net)
		l <- norm_coords(l, ymin=-0.8, ymax=0.9, xmin=-0.9, xmax=0.9)

		plot(net, vertex.frame.color="white", rescale=F, layout=l)

		legend(x = -1.00, y = -0.88, c("No Mutations","Some Mutations", "Frequent Mutations"), pch=21,
	    	col="#777777", pt.cex=c(4, 5, 6), x.intersp=1, y.intersp=1.3, cex=5, bty="n", ncol=1)

		legend(x = -0.40, y = -0.91, c("Non-Linker","Linker"), pch=21,
	    	col="#777777", pt.bg=colrs, pt.cex=5, x.intersp=1.3, y.intersp=1.3, cex=5, bty="n", ncol=1)

		legend(x = 0.1, y = -0.91, c("Non-Kinase","Kinase"), pch=c(21, 22),
	    	col="#777777", pt.cex=5, x.intersp=1, y.intersp=1.3, cex=5, bty="n", ncol=1)

		legend(x = 0.5, y = -0.91, c("Non-Census Gene","Census Gene"), pch=NA,
	    	col="#777777", text.col =label.colrs, pt.cex=5, x.intersp=1.3, y.intersp=1.3, cex=5, bty="n", ncol=1)
	}
}
################################################################################################################################################
formatted.muts <- FormatMutationsbyCancerType(muts, modules)
saveRDS(formatted.muts, formatted.muts.file)

census.genes <- FormatCensusGenes(census)
kinome.genes <- FormatKinomeGenes(kinome)

links <- GetiGraphNetworkLinks(modules, network)
saveRDS(links, links.file)
nodes <- GetiGraphNetworkNodes(nonlinkers.added=modules, network, links, census.genes, kinome.genes, formatted.muts)
saveRDS(nodes, nodes.file)
DrawiGraphModules(links, nodes, modules=module.file)
dev.off()