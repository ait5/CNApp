
Ht.Chroms.RegionProfile <- function(segmented_file) {

  par_chroms <- paste("chr", seq(2, 22, by=2), sep="")
  inpar_chroms <- paste("chr", seq(1, 21, by=2), sep="")
  v_values <- rep(NA, nrow(segmented_file))
  v_colors <- rep(NA, nrow(segmented_file))
  v_values[which(segmented_file[,1]%in%par_chroms)] <- 2
  v_values[which(segmented_file[,1]%in%inpar_chroms)] <- 1
  v_colors[which(segmented_file[,1]%in%par_chroms)] <- "black"
  v_colors[which(segmented_file[,1]%in%inpar_chroms)] <- "gray"


  rr <- segmented_file
  chroms_seq <- paste("chr", seq(1:22), sep="")
  tick_vals_y <- as.vector(do.call(cbind,lapply(1:22, function(x){
    chrom <- chroms_seq[x]
    mat <- rr[which(rr[,1]==chrom),]
    as.character(mat[floor(nrow(mat)/2),4])
  })))
  tick_text_y <- chroms_seq

  y <- rev(as.character(segmented_file[,4]))
  x <- "gplot"
  z <- as.matrix(c(rev(v_values)))
  v_colors <- c(v_colors)

  p_chroms <-  plot_ly(y = y, z = z, x = x, type = "heatmap", colors=v_colors, showscale=F,  text = as.matrix(y), hoverinfo='text')%>%
  layout(
    font = list(size=10),
    xaxis=list(visible=F),
    yaxis=list(visible=T, tickmode="array",
    tickvals=rev(tick_vals_y),
    ticktext=rev(tick_text_y))
  )
  p_chroms
}
