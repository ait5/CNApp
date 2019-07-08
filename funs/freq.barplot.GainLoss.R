

freq.barplot.GainLoss <- function(data, colors, segmented_file, bar_gap=NULL){
  # Frequency plot (as percentages [%])
  # - 'data' has 3 columns:  'loss' ; 'normal' ; 'gain'
  #     'rownames(data)' must be regions (or names) to be x-axis
  # - 'colors' is a vector with many colors as columns there are
  # - 'segmented_file' is a data.frame with chromosome, start, end and region_name information (no header needed)
  # - 'bar_gap' is the space between bars in plot.


  #Preparing new data:
  data$loss <- as.numeric(as.character(data[,"loss"]))
  data$gain <- as.numeric(as.character(data[,"gain"]))
  max_loss <- max(data$loss)
  max_gain <- max(data$gain)
  max_fq <- max(max_loss, max_gain)

  max_term <- 100
  tick_vals <- c(0, 25, 50, 75, 100, 125, 150, 175, 200)
  tick_text <- c("100% -", "75% -", "50% -", "25% -", "0 -", "25% -", "50% -", "75% -", "100% -")

  # if (max_fq > 0 & max_fq < 12.5 ) {
  #   max_term <- 12.5
  #   tick_vals <- c(0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25)
  #   tick_text <- c("12.5%  ", "10%  ", "7.5%  ", "5%  " ,"2.5%  ", "0  ", "2.5%  ", "5%  ",  "7.5%  ", "10%  ", "12.5%  ")
  # } else if (max_fq > 12.5 & max_fq < 25 ) {
  #   max_term <- 25
  #   tick_vals <- c(0, 12.5, 25, 37.5, 50)
  #   tick_text <- c("25%  ",  "12.5%  ", "0  ",  "12.5%  ", "25%  ")
  # } else if (max_fq > 25 & max_fq < 35 ) {
  #   max_term <- 35
  #   tick_vals <- c(0, 17.5, 35, 52.5, 70)
  #   tick_text <- c("35%  ", "17.5%  ", "0  ", "17.5%  ", "35%  ")
  # } else if (max_fq > 35 & max_fq < 50 ) {
  #   max_term <- 50
  #   tick_vals <- c(0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100)
  #   tick_text <- c("50%  ", "37.5%  ", "25%  ", "12.5%  ", "0  ", "12.5%  ", "25%  ", "37.5%  ", "50%  ")
  # } else if (max_fq > 50 & max_fq < 75 ) {
  #   max_term <- 75
  #   tick_vals <- c(0, 25, 50, 75, 100, 125, 150)
  #   tick_text <- c("75%  ", "50%  ", "25%  ", "0  ", "25%  ", "50%  ", "75%  ")
  # } else if (max_fq > 75 & max_fq < 100 ) {
  #   max_term <- 100
  #   tick_vals <- c(0, 25, 50, 75, 100, 125, 150, 175, 200)
  #   tick_text <- c("100%  ", "75%  ", "50%  ", "25%  ", "0  ", "25%  ", "50%  ", "75%  ", "100%  ")
  # } else {
  #   max_term <- 100
  #   tick_vals <- c(0, 25, 50, 75, 100, 125, 150, 175, 200)
  #   tick_text <- c("100%  ", "75%  ", "50%  ", "25%  ", "0  ", "25%  ", "50%  ", "75%  ", "100%  ")
  #
  # }


  data$normal_loss <- max_term - data$loss
  data$normal_gain <- max_term - data$gain


  new_data <- data[,c("normal_loss", "loss", "gain", "normal_gain")]

  n_times <- ceiling(nrow(new_data)/88)
  df <- cbind(new_data, segmented_file[,1])
  n_chroms <- length(unique(segmented_file[,1]))
  l_df <- lapply(1:n_chroms, function(x){
    chrom <- as.character(unique(segmented_file[,1])[x])
    if ( x!=n_chroms ){
      new_row <- c(0, 0, 0, 0, paste(chrom, "_end", sep=""))
      new_mat <- t(replicate(n_times, new_row))
      colnames(new_mat) <- colnames(df)
      rownames(new_mat) <- new_mat[,ncol(new_mat)]
      t <- rbind(df[which(df[,ncol(df)]==chrom),], new_mat)
    } else { t <- df[which(df[,ncol(df)]==chrom),] }
    t
  })
  new_data <- as.data.frame(do.call(rbind, l_df))
  seq_chroms_end <- grep("_end", new_data[,ncol(new_data)])
  chroms_end <- new_data[seq_chroms_end,ncol(new_data)]
  # rownames(new_data)[seq_chroms_end] <- chroms_end
  # new_data <- new_data[,-ncol(new_data)]


  text <- new_data[,-ncol(new_data)]
  text$normal_loss <- rownames(text)
  text$normal_loss[seq_chroms_end] <- ""

  text$normal_gain <- rownames(text)
  text$normal_gain[seq_chroms_end] <- ""

  text$loss <- paste(rownames(text), ": ", new_data$loss, "% loss", sep="")
  text$gain <- paste(rownames(text), ": ", new_data$gain, "% gain", sep="")
  #######

  #Preparing colors:
  colors <- c(colors[2], colors[1], colors[3], colors[2])

  #Preparing chroms plot:
  # chroms <- as.character(unique(segmented_file[,1]))
  par_chroms <- paste("chr", seq(2, 24, by=2), sep="")
  inpar_chroms <- paste("chr", seq(1, 23, by=2), sep="")
  v_values <- rep(NA, nrow(new_data))
  v_colors <- rep(NA, nrow(new_data))
  v_values[which(new_data[,ncol(new_data)]%in%par_chroms)] <- 3
  v_values[which(new_data[,ncol(new_data)]%in%inpar_chroms)] <- 1
  v_values[seq_chroms_end] <- 2
  v_colors[which(new_data[,ncol(new_data)]%in%par_chroms)] <- "black"
  v_colors[which(new_data[,ncol(new_data)]%in%inpar_chroms)] <- "gray"
  v_colors[seq_chroms_end] <- "white"

  rr <- as.data.frame(cbind(rownames(new_data), v_values, as.character(new_data[,ncol(new_data)])))
  chroms_seq <- paste("chr", seq(1:n_chroms), sep="")
  tick_vals_x <- as.vector(do.call(cbind,lapply(1:n_chroms, function(x){
    chrom <- chroms_seq[x]
    mat <- rr[which(rr[,3]==chrom),]
    as.character(mat[floor(nrow(mat)/2),1])
  })))
  tick_text_x <- chroms_seq


  x <- as.character(rownames(new_data))
  y <- "gplot"
  z <- t(as.matrix(c(v_values)))
  v_colors <- unique(c(v_colors))

  p_chroms <-  plot_ly(y = y, z = z, x = x, type = "heatmap", colors=v_colors, showscale=F,  text = as.matrix(x), hoverinfo='x')%>%
  layout(
    font = list(size=10),
    xaxis=list(visible=TRUE, tickmode="array",
    tickvals=tick_vals_x,
    ticktext=tick_text_x, tickfont=list(size=12)),
    yaxis=list(visible=FALSE),
    margin=list(b = 100)
  )
  ########################

  #Preparing FQ plot:
  mt <- new_data[,-ncol(new_data)]

  # x <- as.character(rownames(mt))
  x <- as.character(text[,1])



  # p <- plot_ly(mt, x=x, y=mt[,1],
  #   name=colnames(mt)[1],type = "bar",
  #   marker=list(color=colors[1]),
  #   showlegend=F,
  #   text=text[,1],
  #   hoverinfo='text'
  # )

  # for (i in 2:ncol(mt)){
  #   p <- p %>% add_trace(y=mt[,i],
  #     name=colnames(mt)[i],
  #     marker=list(color=colors[i]),
  #     text=text[,i],
  #     hoverinfo='text'
  #   )
  # }

x <- rownames(mt)
y <- colnames(mt)
dd <- as.data.frame(cbind(x, mt))
colnames(dd)[1] <- "chromosome"

  p <- plot_ly(dd, x = ~chromosome, y = ~normal_loss,
 type = 'bar', name = "", 
 marker=list(color=colors[1], showlegend=FALSE),
      text=text[,1],
      hoverinfo='text', showlegend=FALSE) %>%
  add_trace(y = ~loss, name = 'loss', 
    marker=list(color=colors[2]),
      text=text[,2],
      hoverinfo='text', showlegend=FALSE) %>%
  add_trace(y = ~gain, name = 'gain',
    marker=list(color=colors[3]),
      text=text[,3],
      hoverinfo='text', showlegend=FALSE) %>%
  add_trace(y = ~normal_gain, name = '',
    marker=list(color=colors[4]),
      text=text[,4],
      hoverinfo='text', showlegend=FALSE)



  p <- p %>%   layout(
                yaxis =list(title= "",
                            visible=TRUE, tickmode="array",
                            tickvals=tick_vals,
                            ticktext=tick_text,
                            showgrid=F
  ),
  xaxis = list(visible=T, categoryarray = x, categoryorder = "array"),
  barmode="stack",

  font = list(size=10),
  margin = list(l = 100, r = 80, t = 10, b = 100)
)
# for (i in 2:(length(tick_vals)-1) ){
#   p <- p %>% add_segments(x=x[1], xend=x[length(x)],
#                           y=tick_vals[i], yend=tick_vals[i],
#                           opacity = 0.5,
#                           hoverinfo = "none",
#                           line=list(color="gray", size=1, width=0.5, dash="dash"))
# }
##############

p2 <- subplot(p, p_chroms, nrows=2, shareX=T, heights=c(0.95,0.05))

# When 'bargap' is specified:
if (!is.null(bar_gap)){
   p2$x$layout$bargap <- bar_gap
}

p2
}


