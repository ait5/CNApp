

dendro.Ht.CorrSamples.RegionProfile <- function(data, corr.method, mat_vars, variables_info_mat,
  track_vars, gain.cutoff, loss.cutoff){
    ## App fixed terms ##
    dir_funs <- "funs"
    # Preparing objects needed:
    mat <- data

    vars_dendro <- track_vars
    mat_variables <- mat_vars

    l_vars_annot_in_plot <-  which(colnames(mat_variables)%in%vars_dendro)
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
      #colors_groups <- brewer.pal(length(pca_groups),palette)
      if (var == "BCS") {
        colors_groups <- colorRampPalette(c("deepskyblue", "deepskyblue4"), alpha=TRUE)(length(groups))
      } else if (var == "FCS") {
        colors_groups <- colorRampPalette(c("darkorange", "darkorange4"), alpha=TRUE)(length(groups))
      } else if (var == "GCS") {
        colors_groups <- colorRampPalette(c("darkkhaki", "darkolivegreen"), alpha=TRUE)(length(groups))
      } 
      # else if (var == "ID") {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"PiYG"), alpha=TRUE)(length(groups))
      # } else if (var == "Lineage") {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"Spectral"), alpha=TRUE)(length(groups))
      # } else if ( length(grep("Taylor", var))==1 ) {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"Spectral"), alpha=TRUE)(length(groups))
      # }
      v_colors <- rep(NA, length(groups_vals))
      for (i in 1:length(groups)){
        v_colors[which(groups_vals%in%groups[i])] <- colors_groups[i]
      }

      m_colors <- colors_groups
      names(m_colors) <- groups
      v_empty <- c(v_empty, m_colors)
    }
    #########################

    # Preparing corr samples 'data2':
    cor_samples <- cor(mat, method=corr.method) # 'method_cor' is selected by user (pearson, kendall, spearman)
    cor_samples_mat <- round(cor_samples, 3) # sample correlations matrix

    samples_zero_sd <- which(apply(mat,2,sd)==0)
    if (length(samples_zero_sd)>0){
      data2 <- round(cor(mat[,-samples_zero_sd], method=corr.method), 3)
      annot <- annot[-samples_zero_sd,]

    } else {
      data2 <- round(cor(mat, method=corr.method), 3)
    }
    #########################

    # Preparing 'p_dendro_cor':
    annot <- mat_variables[,l_vars_annot_in_plot]

    if (length(vars_dendro)>1) {
    	rownames(annot) <- as.character(mat_variables[,"ID"])
    }

    p_dendro_cor <- heatmaply(data2, k_col = 1,
      dendrogram="both", plot_method="ggplot",
      col_side_colors=annot,
      col_side_palette=v_empty,
      labCol= NA,
      labRow= NA,
      symm = T,
      hide_colorbar=F,
      margins=c(90, 90, 50, 50))
      #########################

    if (length(track_vars)==1) {
    	p_dendro_cor$x$layout$yaxis2$ticktext[1] <- vars_dendro
    }

      p_dendro_cor
    }
