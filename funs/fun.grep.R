
fun_grep <- function(x, z){
#	HINT: Use it as in 'apply' object-data
#	'x' is a numeric matrix 
#	'y' is a vector of terms to be counted

	e_vector <- rep(NA, length(z))
	for (i in 1:length(z)){
	term <- z[i]
	n_term <- length(which(x==term))
	e_vector[i] <- n_term
	}
	ans <- as.numeric(e_vector)
	names(ans) <- z
	ans
}
     

