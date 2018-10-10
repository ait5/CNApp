

compare.locations.extract.means <- function(i, segmented_file, cna){
	#		'i' = row number of 'segmented_file'

	#		'segmented_file' = first data frame with regions to compare
	#			  chr start     end
	#		1	 chr1     0 2300000

	#		'y' = second data frame with regions to compare
	#			  chr start     end
	#		1	 chr1     0 2300000

	#		'percent' = numeric percentage in order to expand (or shrink) comparision window


	chrom.x <- as.character(segmented_file[i,1])
	start.x <- as.numeric(segmented_file[i,2])
	end.x <- as.numeric(segmented_file[i,3])
	large.x <- end.x-start.x


	# QC to control regions to have positive values in data frame 'y':
	y <- cna
	y.0 <- cna
	y.0[,"qc.positive"] <- as.numeric(y.0[,3]-y.0[,2])
	if (min(y.0[,"qc.positive"])<0){
		rows.bad.value <- which(y.0[,"qc.positive"]<=0)
		y.1 <- y[-rows.bad.value,]

		message.2 <- paste( "Regions ", rows.bad.value, " from 'name df.y' have negative range value and have been ELIMINATED", sep="")
	} else { y.1 <- y }
	###

	if (large.x>0){ #QC-control for regions to have positive values in data frame 'segmented_file'

	#rows in windows


	# start and end point are in window:
	rows.y.1 <- which(y.1[,1]==chrom.x
		&
		( y.1[,2]>=(start.x) & y.1[,2]<=(end.x) & y.1[,3]>=(start.x) & y.1[,3]<=(end.x))
	)

	# start or end point are in window:
	rows.y.2 <- which(y.1[,1]==chrom.x
		&
		( y.1[,2]>=(start.x) & y.1[,2]<=(end.x) | y.1[,3]>=(start.x) & y.1[,3]<=(end.x) )
	)

	rows.y.3 <- which(y.1[,1]==chrom.x
		&
		( y.1[,2]<(start.x) & y.1[,3]>(end.x) )
	)

	rows.y <- c(rows.y.1, rows.y.2, rows.y.3)


	if (length(rows.y)>0){
		#segmented_file[i,"rows.matched.in.y"] <- rows.y
		matched_rows_y <- y.1[rows.y,]

		matched_rows_y$large_in_window <- sapply(1:nrow(matched_rows_y), function(x){
			s.row <- matched_rows_y$loc.start[x]
			e.row <- matched_rows_y$loc.end[x]
			if( (end.x >= s.row & s.row >= start.x) & (end.x >= e.row & e.row >= start.x) ) {
				# start and end are in window
				l <- e.row - s.row
			} else if ( (end.x >= s.row & s.row >= start.x) & (end.x < e.row & e.row >= start.x) ) {
				#start in window
				l <- end.x - s.row
			} else if ( (end.x >= s.row & s.row < start.x) & (end.x >= s.row & s.row >= start.x) ) {
				#end in window
				l <- e.row - start.x
			} else if ( s.row < start.x & e.row > end.x ) {
				#segment embedding window
				l <- large.x
			} else {
				l <- e.row - s.row
			}
			l
		})


		matched_rows_y$percentages <- sapply(1:nrow(matched_rows_y), function(x){
			l <- matched_rows_y$large_in_window[x]
			p <- l/large.x
			if (p>1) {p <- 1}
			p
		}) # percentages from region representing each row

		matched_rows_y$val_with_large <- matched_rows_y$seg.mean * matched_rows_y$percentages

		mean_region <- sum(matched_rows_y$val_with_large)/sum(matched_rows_y$percentages)

	} else if (length(rows.y)==0){
		mean_region <- 0
	}

	mean_region

  } else {
	message.1 <- paste("Region ", i, " from 'name df.x' has negative range value and has NOT being COMPUTED", sep="")
	#		print(message.1)
  }


}
