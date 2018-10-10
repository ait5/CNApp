

genome_build <- "hg38"

NAME <- paste("UCSC_cytoBand_", genome_build, ".txt", sep="")

cytos <- read.table(NAME, sep="\t", header=F)

seq_chroms <- paste("chr", 1:22, sep="")

rows_to_be_selected <- which(cytos[,1]%in%seq_chroms)

sub_cytos <- cytos[rows_to_be_selected,]

order_by_chrom <- order(as.numeric(as.character(gsub("chr", "", sub_cytos[,1]))))

sub_cytos_2 <- sub_cytos[order_by_chrom,]



## Minor cytobands (sub_cytobands)
cyto_name <- "minor_cytobands"
data <- sub_cytos_2[,1:4]
	data[,4] <- paste(data[,1], data[,4], sep="_") 

NAME <- paste("autosomes_", genome_build, "_by_", cyto_name, ".txt", sep="")
write.table(data, NAME, sep="\t", col.names=F, row.names=F, quote=F)



## Major cytobands 
cyto_name <- "major_cytobands"

data <- sub_cytos_2[,1:4]
rownames(data) <- 1:nrow(data)
data$new <- NA
	
for (i in 1:nrow(data)){
	print(paste(i, " from ", nrow(data), sep=""))
	element <- as.character(data[i,4])
	if (length(grep("\\.", element))==1){
		element <- strsplit(as.character(element), "\\.")[[1]][1]
	}
	data[i,"new"] <- element
}
data[,5] <- paste(data[,1], data[,5], sep="_") 
m_bands <- unique(data[,5])

new_data <- as.data.frame(matrix(NA, ncol=4, nrow=length(m_bands)))
for (i in 1:length(m_bands)){
	band <- m_bands[i]
	mat_band <- data[which(data[,5]==band),]
	CHROM <- as.character(unique(mat_band[,1]))
	START <- mat_band[1,2]
	END <- mat_band[nrow(mat_band),3]
	new_row <- c(CHROM, START, END, band)
	new_data[i,] <-	new_row
}
new_data[,2] <- as.numeric(as.character(new_data[,2]))
new_data[,3] <- as.numeric(as.character(new_data[,3]))

sub_cytos_3 <- new_data

data <- sub_cytos_3[,1:4]

NAME <- paste("autosomes_", genome_build, "_by_", cyto_name, ".txt", sep="")
write.table(data, NAME, sep="\t", col.names=F, row.names=F, quote=F)



## Arms 
cyto_name <- "arms"
data <- sub_cytos_3[,1:4]

seq_arms <- c("p", "q")

e_list <- list()
for (i in 1:length(seq_chroms)){
	chrom <- seq_chroms[i]
	mat_chrom <- data[which(data[,1]==chrom),]
	
	for (j in 1:length(seq_arms)){
		arm <- seq_arms[j]
		arm_rows <- grep(arm, mat_chrom[,4])
		mat_chrom_arm <- mat_chrom[arm_rows,]
		START <- mat_chrom_arm[1,2]
		END <- mat_chrom_arm[nrow(mat_chrom_arm),3]
		new_row <- c(chrom, START, END, arm)
		
		e_list[[paste(chrom, arm, sep="_")]] <- new_row
	}
}

new_data <- as.data.frame(do.call(rbind, e_list))
new_data[,2] <- as.numeric(as.character(new_data[,2]))
new_data[,3] <- as.numeric(as.character(new_data[,3]))
new_data[,4] <- rownames(new_data)
rownames(new_data) <- 1:nrow(new_data)

arms_data <- new_data
data <- arms_data

NAME <- paste("autosomes_", genome_build, "_by_", cyto_name, ".txt", sep="")
write.table(data, NAME, sep="\t", col.names=F, row.names=F, quote=F)


## Half-arms
cyto_name <- "half_arms"
data <- arms_data


#by_half_arms

e_list <- list()
for (i in 1:dim(data)[1]){
	row <- data[i, 1:3]
	x.chrom <- as.character(row[1,1])


	row.1.1 <- c(as.character(row[1,1]), row[1,2], round((row[1,3]/2),1))
	row.1.2 <- c(as.character(row[1,1]), row.1.1[3], row[1,3])
	e_list[[i]] <- rbind(row.1.1, row.1.2)
	
}
new_data <- as.data.frame(do.call(rbind, e_list))
	new_data[,4] <- paste(1:dim(new_data)[1], "|", new_data[,1], ":", new_data[,2], "-", new_data[,3], sep="")

new_data[,2] <- as.numeric(as.character(new_data[,2]))
new_data[,3] <- as.numeric(as.character(new_data[,3]))
rownames(new_data) <- 1:nrow(new_data)

half_arms_data <- new_data
data <- half_arms_data

NAME <- paste("autosomes_", genome_build, "_by_", cyto_name, ".txt", sep="")
write.table(data, NAME, sep="\t", col.names=F, row.names=F, quote=F)


#### By (custom) segmentation



seq_bys_num <- c(40000000, 20000000, 10000000, 5000000, 1000000)
seq_bys_names <- c("40Mb", "20Mb", "10Mb", "5Mb", "1Mb")

for (xx in 1:length(seq_bys_num)){

	data <- arms_data

	options("scipen"=100, "digits"=4)
	BY <- seq_bys_num[xx]
	cyto_name <- seq_bys_names[xx]

	print(cyto_name)

	pair_rows <- seq(from=2, to=dim(data)[1], by=2)

	df <- data[pair_rows,1:3]
	df[,2] <- 0


	FUN <- function(x, BY){
	print(BY)
		x.chrom <- as.character(x[1])
		if(length(grep("X", x.chrom))==0 | length(grep("Y", x.chrom))==0){
		x.start <- x[2]
		x.end <- x[3]

		range.1 <- seq(x.start, x.end, by=BY)
		range.1.2 <- range.1[-length(range.1)]
		range.2 <- range.1[2:length(range.1)]

		rep.chrom <- rep(x.chrom, times=length(range.2))

		df.ans <- as.data.frame(as.matrix(cbind(rep.chrom, range.1.2, range.2)))
		df.ans[,2] <- as.numeric(as.character(df.ans[,2]))
		df.ans[,3] <- as.numeric(as.character(df.ans[,3]))
	
		df.ans
		}
	}


	l.ANS <- apply(df, 1, function(x) FUN(x, BY=BY))

	new_data <- do.call(rbind, l.ANS)
	new_data[,4] <- paste(1:dim(new_data)[1], "|", new_data[,1], ":", new_data[,2], "-", new_data[,3], sep="")

	new_data[,2] <- as.numeric(as.character(new_data[,2]))
	new_data[,3] <- as.numeric(as.character(new_data[,3]))
	rownames(new_data) <- 1:nrow(new_data)

	data <- new_data

	NAME <- paste("autosomes_", genome_build, "_by_", cyto_name, ".txt", sep="")
	write.table(data, NAME, sep="\t", col.names=F, row.names=F, quote=F)

}


