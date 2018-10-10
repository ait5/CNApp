
colors.Ht.RegionProfile <- function (data, gain.cutoff, loss.cutoff) {
  # Preparing objects needed:
  mat <- data
  #########################

  # Preparing 'mat_ht':
  MIN_mat <- min(apply(mat,1,min))
  MAX_mat <- max(apply(mat,1,max))

  MIN_MAX <- min(c( (MIN_mat*(-1)) , MAX_mat ) )

  mat_ht <- mat
  cap_value <- MIN_MAX
  for (i in 1:nrow(mat_ht)){
    for (j in 1:ncol(mat_ht)){
      #    mat_ht[i,j] <- (mat_ht[i,j]-mean_tt)/sd_tt
      if (mat_ht[i,j]>=loss.cutoff & mat_ht[i,j]<=gain.cutoff) {mat_ht[i,j] <- 0}
      if (mat_ht[i,j]>=cap_value) {mat_ht[i,j] <- cap_value}
      if (mat_ht[i,j]<= (cap_value*(-1)) ) {mat_ht[i,j] <- (cap_value*(-1))}
    }
  }
  #########################

  # Extracting values from mat_ht:
  e_list <- list()
  for (ii in 1:ncol(mat_ht)){
    values <- mat_ht[,ii]
    e_list[[ii]] <- values
  }
  all_values_mat.0 <- as.vector(do.call(cbind, e_list))
  all_values_mat <- unique(all_values_mat.0)

  tt <- sort(all_values_mat)
  ss <- seq(0,1, length.out=length(tt))
  names(tt) <- ss
  tt2 <- seq(min(tt), max(tt), length.out=length(tt))
  names(tt2) <- ss
  #########################

  color.fun <- function(j, gain_th, loss_th){

    med_gain <- (gain_th * 0.95)
    med_loss <- (loss_th * 0.95)

    items <- length(which(j>=med_gain))
    colfunc_gain <- colorRampPalette(c("lightcoral", "red", "firebrick3"))(items)

    items <- length(which(j<med_gain & j>med_loss))
    colfunc_normal <- colorRampPalette(c("white"))(items)

    items <- length(which(j<=loss_th))
    colfunc_loss <- colorRampPalette(c("royalblue4", "royalblue", "dodgerblue"))(items)

    colors <- c(colfunc_loss, colfunc_normal, colfunc_gain)
    colors
  }
    
  v_colors_ht <- color.fun(tt2, gain.cutoff, loss.cutoff)
  v_colors_ht
}
