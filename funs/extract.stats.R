
extract.stats <- function(x, y, w){
#	'x' data.frame or matrix
#	'y' is column name where to select extracting vector
#	'w' is a second number to extract info from max and mins
n <- which(colnames(x)==y)
z <- as.numeric(as.character(x[,n]))
mean_z <- mean(z)
sd_z <- sd(z)
#density_z <- density(z)
density_z <- 0
top3_min_z <- z[order(z)][1:3]
top3_max_z <- rev(z[order(z)])[1:3]
top3_name_min_z <- x[order(z),w][1:3]
top3_name_max_z <- x[rev(order(z)),w][1:3]

ans <- list("params_fun"= list("x"=x,"y"=y,"w"=w), "mean"=mean_z, "sd"=sd_z, "density_data"=density_z, "top3_min"=top3_min_z, "top3_max"=top3_max_z, "top3_name_min"=top3_name_min_z, "top3_name_max"=top3_name_max_z)
ans
}
