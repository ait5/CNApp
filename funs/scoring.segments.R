#
## Function to assess segmentation and extract values from users segments
#
## INPUTS:
# 'sample_names' = sample names extracted from user input
# 'df' = chromosome-start-end matrix from user input
# 'segmented_file' = chromosome-start-end-value-sample matrix from user input
# 'sample_variable' = sample variable name to be able to select the proper column
#

scoring.segments <- function(x, segmented_file){

	## Fixed inputs
	dir_funs <- "funs"
	##
	chr_control <- grep("chr", segmented_file[,1])
	if (length(chr_control)==0){
		segmented_file[,1] <- paste("chr", segmented_file[,1], sep="")
	}

	name_sample <- as.character(unique(x$ID))

	#print(name_sample)

	if (nrow(x)>0){

		cna <- x[,1:4]

		chr_control <- grep("chr", cna[,1])
		if (length(chr_control)==0){
			cna[,1] <- paste("chr", cna[,1], sep="")
		}


		## compare.locations.R ##

		fun_name <- "compare.locations.extract.means" # function name
		my_fun <- paste(fun_name, ".R", sep="") # file function name
		source_fun <- paste(dir_funs, "/", my_fun, sep="")
		source(source_fun) # sourcing fun into R



		list.means <- lapply(1:nrow(segmented_file), compare.locations.extract.means, segmented_file, cna)




		########################


		df.list.means <- do.call(rbind, list.means)
		if (length(which(is.na(df.list.means)))>0){
			df.list.means[which(is.na(df.list.means))] <- 0
		}

		df.panel <- cbind(segmented_file[,1:3], df.list.means)
		colnames(df.panel) <- c("chr", "start", "end", "mean.reg.sample")
		rownames(df.panel) <- paste(name_sample, 1:dim(df.panel)[1], sep="_")


		df.panel
	}

}
