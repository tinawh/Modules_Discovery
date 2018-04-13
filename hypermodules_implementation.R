# manual hypermodules permutations 

require(data.table); require(tidyr); require(dplyr); require(reshape2); 
require(compare); require(stringr)

# constants 
network.number = 10
bash.script <- "#!/bin/bash\nmodule load java/1.6.0_21\njava -jar /.mounts/labs/reimandlab/private/users/thuang/HyperModules_1.0.2_CMD/HyperModulesCMD-1.0.2.jar -n %s -s %s -c %s -S 0 -t logrank > %s.txt"

# files 
biogrid.net.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv"
net.orig <- fread(biogrid.net.file, header=F)
# mutation.data <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_mutation_data.rds")
# clinical.data <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_fixed_clin.rds")

mut.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/mutation_data/"
clin.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/clinical_data/"

unsplit.clinical <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_fixed_clin.rds")
unsplit.mut <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/unsplit_mutation_data.rds")

# file names 
normal.scripts.folder <- "/.mounts/labs/reimandlab/private/users/thuang/bin_3/original_hypermodules_implementations/"
dir.create(normal.scripts.folder)
normal.results.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/original_hypermodules_implementations/"
dir.create(normal.results.folder)
bin.folder <- "/.mounts/labs/reimandlab/private/users/thuang/bin_3/hypermodules_implementations/"
dir.create(bin.folder)
results.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/hypermodules_implementations/"
dir.create(results.folder)
network.folder <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/04-02-18/biogrid_networks/"
dir.create(network.folder, recursive=T)

GetCancerNames <- function(mut.folder, clin.folder) {
	# get vector of cancer names going into hypermodules 

	setwd(mut.folder)
	muts <- unlist(list.files())
	subbed <- gsub("mut_", "", muts)
	mut.cancers <- gsub("_hugo.csv", "", subbed)

	setwd(clin.folder)
	clins <- unlist(list.files())
	clin.cancers <- gsub("_clin.csv", "", clins)

	tryCatch({identical(mut.cancers, clin.cancers)},
		error = function(e) {
			print("cancers not identical")
			})

	return(mut.cancers)
}

WriteOriginalScripts <- function(cancers, biogrid.net.file, mut.folder, clin.folder, normal.results.folder, normal.scripts.folder) {
	# Write scripts for normal biogrid network

	setwd(normal.scripts.folder)
	lapply(cancers, function(cancer) {
		script <- sprintf(bash.script, biogrid.net.file, paste(mut.folder, "mut_", cancer, "_hugo.csv", sep=""), paste(clin.folder, cancer, "_clin.csv", sep=""), paste(normal.results.folder, cancer, sep=""))
		script.file <- paste(cancer, ".sh", sep="")
		writeLines(script, con=script.file)
		})
}

h_MakeCancerFolder <- function(folder, cancer) {
	dir.path <- paste(folder, cancer, sep="")
	dir.create(dir.path)
}

MakeBinFolders <- function(cancers, h_MakeCancerFolder) {
	# make bin folders per cancer type 

	lapply(cancers, h_MakeCancerFolder, folder=bin.folder)
}

MakeResultsFolders <- function(cancers, h_MakeCancerFolder) {
	# make results folders per cancer type 

	lapply(cancers, h_MakeCancerFolder, folder=results.folder)
}

MakeNetworkFolders <- function(cancers, h_MakeCancerFolder) {
	# make network folders per cancer type 

	lapply(cancers, h_MakeCancerFolder, folder=network.folder)
}

h_GenerateRandomNetwork <- function(network.name, net) {
	# generate random network

	lst <- unique(unlist(net))
	replacement <- data.frame(from=lst, to=sample(lst), stringsAsFactors=F)
	
	random <- net %>% 
		inner_join(replacement, by = c("V1" = "from")) %>% 
		select(to, V2) %>% 
		rename(V1=to) %>% 
		inner_join(replacement, by = c("V2" = "from")) %>% 
		select(to, V1) %>% 
		rename(V2=to)

	# write random network 
	write.table(random, network.name, quote=F, sep="\t", row.names=F, col.names=F)
}

h_WriteRandomNetworksPerCancer <- function(cancer, h_GenerateRandomNetwork, network.folder, network.number, net=net.orig) {
	# generate network.number number of random networks per cancer type 

	seqs <- seq(1:network.number)
	network.names <- paste(paste(paste(network.folder, paste(cancer, "/", sep=""), cancer, sep =""), seqs, sep="_"), ".tsv", sep="")
	lapply(network.names, h_GenerateRandomNetwork, net=net.orig)
}

WriteRandomNetworksForCancers <- function(cancers, network.number=network.number, network.folder=network.folder) {
	lapply(cancers, h_WriteRandomNetworksPerCancer, network.number=network.number, network.folder=network.folder, net.orig)
}

WriteBashScriptsForCancers <- function(cancers) {
	lapply(cancers, h_WriteBashScriptsPerCancer)
}

h_WriteBashScriptsPerCancer <- function(cancer) {
	# Create cancer folder in bin and results then write script per cancer 

	network.folder.cancer <- paste(network.folder, cancer, "/", sep="")
	results.folder.cancer <- paste(results.folder, cancer, "/", sep="")
	bin.folder.cancer <- paste(bin.folder, cancer, "/", sep="")
	setwd(bin.folder.cancer)
	network.lst <- paste(network.folder.cancer, list.files(network.folder.cancer), sep="")

	lapply(network.lst, h_WriteBashScript, results.folder.cancer, bin.folder.cancer)
}

h_WriteBashScript <- function(network.file, results.folder.cancer, bin.folder.cancer) {
	# Write script per network per cancer type 
	netw <- gsub(".tsv", "", unlist(regmatches(network.file, gregexpr("([^/]+$)", network.file, ignore.case=TRUE))))
	cancer.type <- gsub("_[0-9]+", "", netw)
	mut.file <- paste(mut.folder, "mut_", cancer.type, "_hugo.csv", sep="")
	clin.file <- paste(clin.folder, cancer.type, "_clin.csv", sep="") 
	result.file <- paste(results.folder.cancer, netw, sep="")
	script.file <- paste(bin.folder.cancer, netw, ".sh", sep="")
	script.info <- sprintf(bash.script, network.file, mut.file, clin.file, result.file)

	writeLines(script.info, con=script.file)
}

############################################################################################################
cancers <- GetCancerNames(mut.folder, clin.folder)
WriteOriginalScripts(cancers, biogrid.net.file, mut.folder, clin.folder, normal.results.folder, normal.scripts.folder)
MakeBinFolders(cancers, h_MakeCancerFolder)
MakeResultsFolders(cancers, h_MakeCancerFolder)
MakeNetworkFolders(cancers, h_MakeCancerFolder)

# network.lst <- MakeNetworklist(network.number, loc)
WriteBashScriptsForCancers(cancers)

