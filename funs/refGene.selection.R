#
## Function to select segmentation file
#
## INPUTS:
# 'segmentation_name' = segmentation type selected by user to be used as genome segmentation
# 'genome_build' = genome build version (options: 'hg19' (build 37) or 'hg38')
#
#

refGene.selection <- function(genome_build){

	aux_files_dir <- "aux_files"
	name_dir_hg <- paste("genes_in_", genome_build, sep="")
	name_file_hg <- paste(aux_files_dir, "/", name_dir_hg, "/", paste("refGene_", genome_build, ".txt", sep=""), sep="")
	ans <- read.table(name_file_hg, sep="\t", header=T)
}

