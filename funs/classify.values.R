
classify.values <- function(x, gain, loss){
#	HINT: Use 'evaluate.region.values' as 'lapply' object-data
#	'x' is a numeric term from a list
#	'y' is a threshold value to classify list values to UP, DOWN or EQUAL

	if (!is.na(x)){
		if (x <= loss){
			ans <- "loss"
		} else if (x >= gain){
#			if (x >= high.gain){
#				ans <- "high-gain"
#			} else {
				ans <- "gain"
#			}	
		} else {
			ans <- "normal"
		}
	} else {ans <- NA}
	ans
}
