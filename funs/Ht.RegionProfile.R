

Ht.RegionProfile <- function (data, mat_vars, segmented_file, variables_info_mat,
  track_vars, gain.cutoff, loss.cutoff) {
    ## App fixed terms ##
    dir_funs <- "funs"
    # Preparing objects needed:
    mat <- data
    colnames(mat) <- paste(" -", colnames(mat), "- ", sep="")
    sample_names <- colnames(mat)
    seg_regions_names <- rownames(mat)

    var_tracks <- track_vars
    mat_variables <- mat_vars
    mat_annot <- as.matrix(t(mat_variables))

    re_order_seq <- seq(nrow(mat),1, length=nrow(mat))

    l_vars_annot_in_plot <-  which(colnames(mat_variables)%in%var_tracks)
    #########################

    # Preparing heatmap colors:
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

    # Preparing hetmap colors:
    fun_name <- "colors.Ht.RegionProfile" # function name
    my_fun <- paste(fun_name, ".R", sep="") # file function name
    source_fun <- paste(dir_funs, "/", my_fun, sep="")
    source(source_fun) # sourcing fun into R

    v_colors_ht <- colors.Ht.RegionProfile(data, gain.cutoff, loss.cutoff)
    ############################################



    p_cna_profile <- plot_ly(y=rownames(mat[re_order_seq,]),
    x=colnames(mat),
    z = as.matrix(mat_ht[re_order_seq,]),
    type = "heatmap",
    colors=v_colors_ht,
    text = round(as.matrix(mat[re_order_seq,]),3),
    hoverinfo='x+y+text')%>%
    colorbar(
      thicknessmode="fraction",
      thickness=0.02,


      tickmode="array",
      tickvals=c(((cap_value*(-1))*0.8), cap_value*0.8),
      ticktext=c("loss", "gain"),
      ticklen=0
    )%>%
    layout(
      font = list(size=10),
      margin = list(l = 80, r = 80, t = 40, b = 40),
      xaxis=list(visible=F),
      yaxis=list(visible=F)
    )
    p_cna_profile

  }
