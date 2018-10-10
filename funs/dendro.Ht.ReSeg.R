

dendro.Ht.ReSeg <- function(data, dendro, mat_vars, variables_info_mat, fontsize_col=NULL){
    ## App fixed terms ##
    dir_funs <- "funs"
    # Preparing objects needed:
    mat <- data

    # vars_dendro <- track_vars
    mat_variables <- mat_vars

    # l_vars_annot_in_plot <-  which(colnames(mat_variables)%in%vars_dendro)
    if (is.null(fontsize_col)){
      fontsize_col <- 10
    }
    #########################

    # Preparing colors for annotation tracks:
    v_empty <- vector()
    for (j in 1:ncol(mat_variables)){
      var <- colnames(mat_variables)[j]
      groups_vals <- mat_variables[,which(colnames(mat_variables)==var)]
      groups <- sort(unique(groups_vals))
      class_var <- variables_info_mat[which(variables_info_mat[,1]==var),"class_var"]
      palette <- as.character(variables_info_mat[which(variables_info_mat[,1]==var),"color_palette"])

      colors_groups <- colorRampPalette(brewer.pal(11,palette))(length(groups))
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



    # Preparing 'p_dendro':

    annot <- mat_variables
    rownames(annot) <- as.character(mat_variables[,"ID"])

    ht <- heatmaply(t(mat), k_col = 1, k_row=1,
      dendrogram=dendro, plot_method="ggplot",
      row_side_colors=annot,
      row_side_palette=v_empty,
      hide_colorbar=T,
      fontsize_col= fontsize_col,
      labRow= NA,
      margins=c(150, 100, 100, 100),
    cexCol=NA,
  subplot_widths=c(0.4,0.35,0.25))

      ht
    }
