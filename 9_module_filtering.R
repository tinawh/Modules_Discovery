# filter hypermodules results

require(data.table); require(tidyr); require(dplyr); require(reshape2); require(compare)
require(survival); require(survminer); require(tcR); require(purrr); require(gplots)

options(warn = -1)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-05-18/"
dir.create(loc)

# constants
hypermodules.dir <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/hypermodules_biogrid_10_1/"
hypermodules.prefix <- "hypermodules_10_" # number of permutations
patient.threshold <- 10 # threshold for number of patients in a module
coxph.threshold <- 0.05
jaccard.cutoff <- 0.2


# files
clins <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/unsplit_fixed_clin.rds")

# file names
nosubset.results.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-05-18/nosubset_hypermodules_results_permut_10.rds"
heatmap.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_3/03-05-18/module_clustering_heatmap.pdf"

ReadHypermodulesResults <- function(hypermodules.dir, hypermodules.prefix) {
	# read in hypermodules results and remove cancers without results

	setwd(hypermodules.dir)
	file.paths <- list.files()

	lst <- lapply(file.paths, function(fil) {
		if (!file.size(fil) == 0) {
			fread(fil, skip = 3, verbose = TRUE)
			}
		})

	nam <- gsub(hypermodules.prefix, "", file.paths)
	names(lst) <- (gsub(".txt", "", nam)) 

	# remove NULLs
	lst[sapply(lst, is.null)] <- NULL

	for (i in names(lst)) {
		if (ncol(lst[[i]]) == 2) {
			lst[[i]] <- NULL
		}
	}

	binded <- bind_rows(lst, .id = "cancer")

	return(binded)
}

RemoveSmallModules <- function(hypermodules.results, patient.threshold) {
	# remove modules that have less than patient.threshold patients
	results <- hypermodules.results %>%
					filter(Number_patients >= patient.threshold) 
	return(results)
}

SplitHypermodulesResults <- function(hypermodules.results) {
	results.lst <- split(hypermodules.results, hypermodules.results$cancer)
	return(results.lst)
}

SplitClinicalResults <- function(clins, cancer.names = names(hypermodules.results.lst)) {
	clins.ab <- clins %>% 
					filter(Study_Name %in% cancer.names)

	splitted <- split(clins.ab, clins.ab$Study_Name)
	return(splitted)
}

h_CalculateLogRankCoxph <- function(result, clin) {
	# returns logrank and coxph for each module in a cancer type 
	splitted <- split(result, result$Module)
	split.rows <- list()
	for (df in names(splitted)) {
		df_update <- splitted[[df]] %>%
						mutate(List_of_patients = strsplit(List_of_patients, split = ",")) %>% 
						select(Seed, Module, List_of_patients) %>%
						unnest(List_of_patients) %>%
						unique() %>%
						mutate(is_module_mutated = 1) %>%
						full_join(clin, by = c("List_of_patients" = "patient")) %>%
						mutate(vital_status = gsub("DECEASED", 1, vital_status)) %>%
						mutate(vital_status = gsub("ALIVE", 0, vital_status)) %>%
						mutate(vital_status = as.numeric(vital_status)) %>% 
						select(List_of_patients, is_module_mutated, vital_status, surv_time)

					df_update$is_module_mutated[is.na(df_update$is_module_mutated)] <- 0
					split.rows[[df]] <- df_update
	}
	# calculate logrank 
	lst = list()
	for (df in names(split.rows)) {
		tryCatch(sdiff <- survdiff(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]]), warning = function(w) print(paste(df, "sdiff error")))
		sdiff <- survdiff(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]])
		p_logrank = 1-pchisq(sdiff$chisq, df=1)
		lst[[df]] <- p_logrank
	}

	chr <- as.numeric(as.character(lst))
	names(chr) <- names(lst)

	# calculate coxph & survival correlation
	lst_c <- list()
	lst_coeff <- list() # correlation
	for (df in names(split.rows)) {
		tryCatch(h0 <- coxph(Surv(surv_time, vital_status) ~ 1, data = split.rows[[df]]), warning = function(w) print(paste(df, "h0 error")))
		h0 <- coxph(Surv(surv_time, vital_status) ~ 1, data = split.rows[[df]])
		tryCatch(h1 <- coxph(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]]), warning = function(w) print(paste(df, "h1 error")))
		h1 <- coxph(Surv(surv_time, vital_status) ~ is_module_mutated, data = split.rows[[df]])
		p_cox <- anova(h0, h1)[2,4] 
		lst_c[[df]] <- p_cox

		# correlation
		coeff <- h1[[1]]
		names(coeff) <- NULL
		lst_coeff[[df]] <- coeff 
	}

	chr_c <- as.numeric(as.character(lst_c))
	names(chr_c) <- names(lst_c)

	chr_c_corr <- as.numeric(as.character(lst_coeff))
	names(chr_c_corr) <- names(lst_coeff)

	dat <- data.frame(chr, chr_c, chr_c_corr) %>%
				tibble::rownames_to_column() %>%
				rename(module = rowname, logrank = chr, coxph = chr_c, coxph_corr = chr_c_corr) %>%
				inner_join(result, by = c("module" = "Module"))
	return(dat)
}


RemoveNonsignificantCoxph <- function(test.results, coxph.threshold) {
	# Remove modules with coxph greater than 0.05

	filtered <- bind_rows(test.results) %>% 
					filter(coxph <= coxph.threshold) 
return(filtered)		
}

CreateSubsetMatrixPerPatientList <- function(n, h.results) {
	new <- h.results %>% 
		mutate(List_of_patients = strsplit(List_of_patients, split=","))
	patients <- new$List_of_patients
	lapply(patients, function(patient.list) {
		all(patient.list %in% patients[[n]])
		})
}

CreateMatrixPerCancer <- function(coxphfiltered.results) {
	splitted <- split(coxphfiltered.results, coxphfiltered.results$cancer)
	cancers <- list()
	for (canc in names(splitted)) {
		cancers[[canc]] <- do.call(rbind, lapply(seq(nrow(splitted[[canc]])), CreateSubsetMatrixPerPatientList, splitted[[canc]]))
	}
	return(cancers)
}

FilterMatrix <- function(cancer.matrices, splitted) {
	new.lst <- list()
	for (canc in names(cancer.matrices)) {
			matx <- cancer.matrices[[canc]]
			num.modules <- seq(dim(matx)[1])
			patient.nums <- splitted[[canc]]$Number_patients 
			for(i in num.modules) {matx[i,i] <- FALSE}
			lst <- apply(which(matx==TRUE, arr.ind=TRUE), 1, as.list)
			dups <- sapply(lst, function(sett) { 
				ifelse(patient.nums[sett[[1]]] > patient.nums[sett[[2]]], sett[[2]], sett[[1]])
			})
		new.lst[[canc]] <- splitted[[canc]][setdiff(num.modules, dups), ]
	}
	return(bind_rows(new.lst))
}

# RemoveRedundantPatientModules <- function(coxphfiltered.results) {
# 	# Remove modules that are subsets of larger modules

# 	names(coxphfiltered.results)[1] <- "Module" # rename
# 	coxphfiltered.results <- split(coxphfiltered.results, f = coxphfiltered.results$cancer)

# 	filtered.results <- list()

# 	for (result in names(coxphfiltered.results)) {
# 		splitted <- coxphfiltered.results[[result]] %>% 
# 						mutate(List_of_patients = strsplit(List_of_patients, split = ","))

# 		patient.list <- splitted$List_of_patients
# 		names(patient.list) <- coxphfiltered.results[[result]]$Module

# 		sorted <- character()
# 		for (lst in names(patient.list)) {
# 			ss <- sort(patient.list[[lst]])
# 			sorted[[lst]] <- paste(ss, collapse = ",")
# 		}

# 		splitted %>% 
# 			mutate(List_of_patients = sorted) %>% 
# 			mutate(List_of_patients = strsplit(List_of_patients, split = ",")) -> sorted.splitted

# 		sorted.test <- sorted.splitted$List_of_patients
# 		names(sorted.test) <- sorted.splitted$Module

# 		uniques <- unique(sorted.test)
# 			if (compare::compare(sorted.test, uniques, allowAll = TRUE)[1] == FALSE) {
# 					dups <- sorted.test[names(sorted.test[which(sorted.test[duplicated(sorted.test)] %in% sorted.test)])]
# 				} else {
# 					dups <- list()
# 				}
			
# 		sortings <- sorted.test
# 		for (entry in names(sorted.test)) {
# 			 x <- sorted.test[[entry]]
# 			 deleted <- sortings
# 			 deleted[[entry]] <- NULL
# 			 for (y in names(deleted)) {
# 			 	if (all(x %in% deleted[[y]])) {
# 			 		sortings[[y]] <- NULL
# 			 	}
# 			}
# 		}

# 		new.lst <- c(sortings, dups)

# 		screened.list <- list()
# 		for (mod in names(new.lst)) {
# 		screened.list[[mod]] <- paste(new.lst[[mod]], collapse = ",")
# 		}

# 		new.dframe <- data.frame(Module = as.character(names(screened.list)), stringsAsFactors = F) %>% 
# 						inner_join(coxphfiltered.results[[result]], by = "Module")

# 		new.dframe <- new.dframe[, c("Module", "logrank", "coxph", "cancer", "coxph_corr", "Number_patients", "List_of_patients")]

# 		filtered.results[[result]] <- new.dframe
# 	}
# 	binded.results <- bind_rows(filtered.results)

# 	return(binded.results)	
# }

Splitting <- function(removedsubset.results) {
	splitted <- split(removedsubset.results, f=removedsubset.results$cancer) 

	return(splitted)
}

h_RemoveSimilarModules <- function(n, orig.splitted.cancer) {

	 lst <- as.list(orig.splitted.cancer$List_of_patients)
	 coeff <- lapply(lst, function(x) {
	 	jaccard.index(unlist(strsplit(x, split = ",")), unlist(strsplit(lst[[n]], split = ",")))
	 	})

	 return(unlist(coeff))
}

GetMatrix <- function(splitted, h_RemoveSimilarModules) {

	split.lst <- list() 

	for (cancer in names(splitted)) {
		named.lst <- lapply(1:nrow(splitted[[cancer]]), h_RemoveSimilarModules, splitted[[cancer]])
		# names(named.lst) <- splitted[[cancer]]$Module
		split.lst[[cancer]] <- do.call(rbind, named.lst)
	}

	return(split.lst)
}

FilterMatrix <- function(split.lst.cancer, orig.splitted.cancer, jaccard.cutoff) {

	sample_size <- orig.splitted.cancer$Number_patients 
	sig_factor <- orig.splitted.cancer$coxph 
	modules <- 1:length(sample_size) 

	retrieve_less_sig_set = function(i, j){
	  if(sample_size[[i]] == sample_size[[j]]){
	    return(ifelse(sig_factor[[i]] > sig_factor[[j]],j,i))
	  }else if(sample_size[[i]] > sample_size[[j]])
	    return(j)
	  else(return(i))
	}

	for(i in 1:length(sample_size)){
	  for(j in 1:length(sample_size)){
	    if(split.lst.cancer[i,j]  > jaccard.cutoff & i != j){
	      modules <- setdiff(modules, retrieve_less_sig_set(i, j))
	    }
	  }
	}
	return(data.frame(orig.splitted.cancer[modules,], stringsAsFactors=F))
}


GetDistanceMatrix <- function(splitted, h_RemoveSimilarModules) {

	split.lst <- list() 

	for (cancer in names(splitted)) {
		named.lst <- lapply(1:nrow(splitted[[cancer]]), h_RemoveSimilarModules, splitted[[cancer]])
		# names(named.lst) <- splitted[[cancer]]$Module
		split.lst[[cancer]] <- 1 - do.call(rbind, named.lst)
	}

	return(split.lst)
}

h_GetFET <- function(n, orig.splitted.cancer) {

	 lst <- as.list(orig.splitted.cancer$List_of_patients)
	 names(lst) <- orig.splitted.cancer$Module
	 new.lst <- list()
	 for (x in names(lst)) {
	 	all.patients <- union(unlist(strsplit(lst[[x]], split = ",")), unlist(strsplit(lst[[n]], split = ",")))
	 	if (all(unlist(strsplit(lst[[x]], split = ",")) == unlist(strsplit(lst[[n]], split=",")))) {
	 		1
	 	} else {
	 		tryCatch(
	 		{new.lst[[x]] <- fisher.test(all.patients %in% unlist(strsplit(lst[[x]], split = ",")), all.patients %in% unlist(strsplit(lst[[n]], split = ",")))$p.value
	 		}, error = function(e) {print(paste(unique(orig.splitted.cancer$cancer), lst[[x]], n))}
	 		)
	 	} 
	 }

	#  coeff <- lapply(lst, function(x) {
	#  	all.patients <- union(unlist(strsplit(x, split = ",")), unlist(strsplit(lst[[n]], split = ",")))
	#  	if (all(unlist(strsplit(x, split = ",")) == unlist(strsplit(lst[[n]], split=",")))) {
	#  		0
	#  	} else {
	#  		tryCatch(
	#  		{fisher.test(all.patients %in% unlist(strsplit(x, split = ",")), all.patients %in% unlist(strsplit(lst[[n]], split = ",")))$p.value
	#  		}, error = function(e) {print(paste(unique(orig.splitted.cancer$cancer), x, n))}
	#  		)
	#  	# jaccard.index(unlist(strsplit(x, split = ",")), unlist(strsplit(lst[[n]], split = ",")))
	#  }
	# })

	 return(new.lst)
}

GetFETMatrix <- function(splitted, h_GetFET) {

	split.lst <- list() 

	for (cancer in names(splitted)) {
		named.lst <- lapply(1:nrow(splitted[[cancer]]), h_GetFET, splitted[[cancer]])
		# names(named.lst) <- splitted[[cancer]]$Module
		split.lst[[cancer]] <- do.call(rbind, named.lst)
	}

	return(split.lst)
}

DrawHeatMap <- function(distance.lst, heatmap.file) {
	
	pdf(heatmap.file)
	distance.lst.new <- distance.lst[which(lapply(distance.lst, length) > 1)]

	for (cancer.m in names(distance.lst.new)) {
		rc <- rainbow(nrow(distance.lst.new[[cancer.m]]), start=0, end=.3)
 		cc <- rainbow(ncol(distance.lst.new[[cancer.m]]), start=0, end=.3)
 		heatmap.2(distance.lst.new[[cancer.m]], trace="none", main=cancer.m) 
	}
	# lapply(distance.lst.new, function(cancer.m) {
	# 	rc <- rainbow(nrow(cancer.m), start=0, end=.3)
 # 		cc <- rainbow(ncol(cancer.m), start=0, end=.3)
 # 		heatmap.2(cancer.m)
	# 	})
}

GetSurvivalCorrelations <- function(binded) {
	# Get number of positive and negative survival correlations

	wt.dobetter <- sum(binded$coxph_corr > 0)
	mut.dobetter <- sum(binded$coxph_corr < 0)

	nums <- list(wt.dobetter, mut.dobetter)
	names(nums) <- c("wt_dobetter", "mut_dobetter")

	return(nums)
}

########################################################################################################
hypermodules.results <- ReadHypermodulesResults(hypermodules.dir, hypermodules.prefix)
removed.results <- RemoveSmallModules(hypermodules.results, patient.threshold)

hypermodules.results.lst <- SplitHypermodulesResults(removed.results)
clins.lst <- SplitClinicalResults(clins)

test.results <- list()
	for (i in names(hypermodules.results.lst)) {
		result <- h_CalculateLogRankCoxph(hypermodules.results.lst[[i]], clins.lst[[i]])
		test.results[[i]] <- result
	}

coxphfiltered.results <- RemoveNonsignificantCoxph(test.results, coxph.threshold)

cancer.matrices <- CreateMatrixPerCancer(coxphfiltered.results)
removedsubset.results <- FilterMatrix(cancer.matrices, splitted)

# removedsubset.results <- RemoveRedundantPatientModules(coxphfiltered.results)

splitted <- Splitting(removedsubset.results)

split.lst <- GetMatrix(splitted, h_RemoveSimilarModules)

distance.lst <- GetDistanceMatrix(splitted, h_GetJaccardDistance) 

FET.m <- GetFETMatrix(splitted, h_GetFET) 

DrawHeatMap(distance.lst, heatmap.file)
dev.off()

nosubset.results  <- list()
for (cancer in names(splitted)) {
	nosubset.results[[cancer]] <- FilterMatrix(split.lst[[cancer]], splitted[[cancer]], jaccard.cutoff)
}

binded <- bind_rows(nosubset.results) # bind to single dataframe 
saveRDS(binded, nosubset.results.file)

survival.correlations <- GetSurvivalCorrelations(binded)



