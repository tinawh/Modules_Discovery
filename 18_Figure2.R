# Figure 2. Sub-networks characterize patient survival in XX cancer types
# show all networks as a cytoscape visualisation. group them spatially into areas for different cancer types
# count modules: have one stacked barplot with two colors, green/orange, positive/negative survival assoc. for each cancer type

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(compare)
require(stringr); require(igraph); require(ggpubr); require(rvest); require(httr); require(RJSONIO)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-19-18/Figure_2/"
dir.create(loc)

# files
modules <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/modules_added_nonlinkers_permut_10.rds")
nodes <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/nodes.rds")
links <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_3/03-12-18/links.rds")

# constants
url <- "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations"
tcga.table <- '//*[@id="node-677"]/div/div[1]/div/div/table'

# file names 
module.distribution.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-19-18/Figure_2/module_distribution.pdf"

GetCancerAbbrev <- function(url, tcga.table) {
	# Get tcga cancer abbreviations

	ids <- url %>% 
		read_html() %>% 
		html_nodes(xpath=tcga.table) %>% 
		html_table() %>% 
		as.data.frame(stringsAsFactors=F) %>% 
		mutate(Study.Name=gsub(" ", "_", Study.Name))

	return(ids)
}

SplitNetworks <- function(nodes, links) {

}

DrawGroupedNetworks <- function() {

	brain.nodes <- nodes[2:27]
	brain.links <- links[2:27]

	pdf(width = 1000, height = 1000)
	par(mfrow=c(4,7), family="sans")
	for (mod in names(brain.nodes)) {
		net <- graph_from_data_frame(d=brain.links[[mod]], vertices=brain.nodes[[mod]], directed=F)
			
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

			# par(family = "sans")

			l <- layout_with_kk(net)
			l <- norm_coords(l, ymin=-0.8, ymax=0.9, xmin=-0.9, xmax=0.9)

			plot(net, vertex.frame.color="white", rescale=F, layout=l)
	}
		legend(x = -1.00, y = -0.88, c("No Mutations","Some Mutations", "Frequent Mutations"), pch=21,
	    	col="#777777", pt.cex=c(4, 5, 6), x.intersp=1, y.intersp=1.3, cex=5, bty="n", ncol=1)

		legend(x = -0.40, y = -0.91, c("Non-Linker","Linker"), pch=21,
	    	col="#777777", pt.bg=colrs, pt.cex=5, x.intersp=1.3, y.intersp=1.3, cex=5, bty="n", ncol=1)

		legend(x = 0.1, y = -0.91, c("Non-Kinase","Kinase"), pch=c(21, 22),
	    	col="#777777", pt.cex=5, x.intersp=1, y.intersp=1.3, cex=5, bty="n", ncol=1)

		legend(x = 0.5, y = -0.91, c("Non-Census Gene","Census Gene"), pch=NA,
	    	col="#777777", text.col =label.colrs, pt.cex=5, x.intersp=1.3, y.intersp=1.3, cex=5, bty="n", ncol=1)
}

DrawStackedBarplot <- function(modules, ids, module.distribution.file) {
	# Draw stacked barplot by cancer type to show survival correlation

	modules$coxph_corr <- ifelse(modules$coxph_corr > 0, "Negative", "Positive")

	selected <- modules %>% 
		inner_join(ids, by=c("cancer" = "Study.Name")) %>%
		ungroup() %>% 
		select(coxph_corr, Study.Abbreviation)

	dframe <- as.data.frame(table(selected), stringsAsFactors=F)
	# dframe <- dframe[order(dframe$coxph_corr, decreasing = T),]

	ordered <- names(sort(table(selected$Study.Abbreviation), decreasing=T))
	caption <- paste("n =", nrow(modules))

	pdf(module.distribution.file)

	ggbarplot(dframe, 
		x="Study.Abbreviation", 
		y="Freq", 
		color="coxph_corr", 
		fill="coxph_corr",  
		palette=c("#a5bc6a", "#fea002"),
		width=0.6, 
		# alpha=0.8,
		legend="bottom",
		legend.title="Correlation",
		# legend.labels=c("Negative", "Positive"),
		order=ordered,
		x.text.angle=45,
		caption=caption) +
	labs(title="Subnetwork Distribution Across Cancers", 
		x="Cancer",
		y="Number of Modules")
}

########################################################################################################################
ids <- GetCancerAbbrev(url, tcga.table)
DrawStackedBarplot(modules, ids, module.distribution.file)
