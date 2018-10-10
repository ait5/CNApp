#
## Function to select segmentation file
#
## INPUTS:
# 'segmentation_name' = segmentation type selected by user to be used as genome segmentation
# 'genome_build' = genome build version (options: 'hg19' (build 37) or 'hg38')
#
#

segmentation.selection <- function(segmentation_name, genome_build){

	aux_files_dir <- "aux_files"
	name_seg_hg <- paste("segmented_files_", genome_build, sep="")
	name_file <- paste(aux_files_dir, "/", name_seg_hg, ".txt", sep="")
	df_segmented_files <- read.table(name_file, sep="\t", header=T)
	
	row <- df_segmented_files[which(df_segmented_files[,"name"]==segmentation_name),]

	in_SEG_FILE <- paste(aux_files_dir, "/", name_seg_hg, "/", row[,1], sep="")
	segmented_file <- read.table(in_SEG_FILE, sep="\t", header=F)
	colnames(segmented_file)[1:4] <- c("chr", "start", "end", "reg_name")
	
	ans <- segmented_file
}

