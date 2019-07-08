

dendro.Ht.RegionProfile <- function(data, dendro, mat_vars, segmented_file, variables_info_mat,
  track_vars, gain.cutoff, loss.cutoff, fontsize_row=NULL){
    ## App fixed terms ##
    dir_funs <- "funs"
    # Preparing objects needed:
    mat <- data

    vars_dendro <- track_vars
    mat_variables <- mat_vars

    l_vars_annot_in_plot <-  which(colnames(mat_variables)%in%vars_dendro)
    if (is.null(fontsize_row)){
      fontsize_row <- 10
    }
    #########################

    # Preparing colors for annotation tracks:
    v_empty <- vector()
    for (j in 1:length(vars_dendro)){
      var <- vars_dendro[j]
      groups_vals <- mat_variables[,which(colnames(mat_variables)==var)]
      groups <- sort(unique(groups_vals))
      class_var <- variables_info_mat[which(variables_info_mat[,1]==var),"class_var"]
      palette <- variables_info_mat[which(variables_info_mat[,1]==var),"color_palette"]

      colors_groups <- colorRampPalette(brewer.pal(11,palette))(length(groups))
      
      if (var == "BCS") {
        colors_groups <- colorRampPalette(c("deepskyblue", "deepskyblue4"), alpha=TRUE)(length(groups))
      } else if (var == "FCS") {
        colors_groups <- colorRampPalette(c("darkorange", "darkorange4"), alpha=TRUE)(length(groups))
      } else if (var == "GCS") {
        colors_groups <-colorRampPalette(c("darkkhaki", "darkolivegreen"), alpha=TRUE)(length(groups))
      } 
      # else if (var == "ID") {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"PiYG"), alpha=TRUE)(length(groups))
      # } else if (var == "Lineage") {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"Spectral"), alpha=TRUE)(length(groups))
      # } else if ( length(grep("Taylor", var))==1 ) {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"Spectral"), alpha=TRUE)(length(groups))
      # }
      
      #colors_groups <- brewer.pal(length(pca_groups),palette)
      v_colors <- rep(NA, length(groups_vals))
      for (i in 1:length(groups)){
        v_colors[which(groups_vals%in%groups[i])] <- colors_groups[i]
      }
      
      m_colors <- colors_groups
      names(m_colors) <- groups
      v_empty <- c(v_empty, m_colors)
    }
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
        # if (mat_ht[i,j]>=loss.cutoff & mat_ht[i,j]<=gain.cutoff) {mat_ht[i,j] <- 0}
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

    # Preparing 'p_dendro':
    rownames(mat_ht) <- as.character(segmented_file[,4])

    annot <- as.matrix(mat_variables[,l_vars_annot_in_plot])
    colnames(annot) <- track_vars
    rownames(annot) <- as.character(mat_variables[,"ID"])

    # if (length(vars_dendro)>1) {
    # 	rownames(annot) <- as.character(mat_variables[,"ID"])
    # }
    

    p_dendro <- heatmaply(mat_ht, k_col = 1,
      dendrogram=dendro, plot_method="ggplot",
      col_side_colors=annot,
      col_side_palette=v_empty,
      #custom_hovertext=text_mat,
      colors=v_colors_ht,
      hide_colorbar=F,
       fontsize_row= fontsize_row,
      labCol= NULL,
      margins=c(40, 80, 40, 80))

    if (length(track_vars)==1) {
    	p_dendro$x$layout$yaxis2$ticktext[1] <- vars_dendro
    }

      p_dendro
    }
