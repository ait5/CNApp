
Ht.Annotation.2.RegionProfile <- function(data, mat_vars,
  variables_info_mat, track_vars) {

    # Preparing objects needed:
    mat <- data
    sample_names <- colnames(mat)
    seg_regions_names <- rownames(mat)

    var_tracks <- track_vars
    mat_variables <- mat_vars
    mat_annot <- as.matrix(t(mat_variables))

    l_vars_annot_in_plot <-  which(colnames(mat_variables)%in%var_tracks)
    seq_1 <- seq(1,length(l_vars_annot_in_plot), by=1)
    #########################

    plotList2 <- function(x){
      i <- l_vars_annot_in_plot[x]
      by_for_colorbar <- 0.05
      y <- colnames(mat_variables)[i]
      class_var <- variables_info_mat[which(variables_info_mat[,1]==y),2]
      palette <- variables_info_mat[which(variables_info_mat[,1]==y),3]
      x <- colnames(mat)

      text <- as.matrix(as.character(mat_annot[i,]))

      if (class_var=="numeric") {
        mat_test <- as.matrix(as.numeric(as.character(mat_annot[i,])))
        tickvals_0 <- unique(as.numeric(as.character(mat_test)))
        tickvals <- c(min(tickvals_0), mean(tickvals_0), max(tickvals_0))
        range <- tickvals[length(tickvals)]-tickvals[1]
        tickvals[1] <- tickvals[1] + (range*0.1)
        tickvals[length(tickvals)] <- tickvals[length(tickvals)] - (range*0.1)
        ticktext_0 <- unique(as.numeric(as.character(text)))
        ticktext <- round(c(min(ticktext_0), mean(ticktext_0), max(ticktext_0)),2)
      } else if (class_var=="categoric") {
        mat_test <- as.matrix(as.numeric(as.factor(mat_annot[i,])))
        tickvals <- unique(as.numeric(as.character(mat_test)))
        range <- tickvals[length(tickvals)]-tickvals[1]
        tickvals[1] <- tickvals[1] + (range*0.1)
        tickvals[length(tickvals)] <- tickvals[length(tickvals)] - (range*0.1)
        ticktext <- unique(as.character(text))
      }

      y <- y
      x <- rep("gplot", 1)
      z <- as.matrix(rep(1, 1))
      v_colors <- "white"
      p2 <- plot_ly(y = y, z = z, x = x, type = "heatmap", colors=v_colors, showscale=F, hoverinfo='none')%>%
      layout(
        font = list(size=10),
        xaxis=list(visible=F),
        yaxis=list(visible=F)
      )
      p2

    }
    if (length(l_vars_annot_in_plot)>1) {
      jj2 <- lapply(seq_1, plotList2)
      p_annotation2 <- subplot(jj2, nrows=length(jj2), shareX=T)
    } else {
      p_annotation2 <- plotList2(seq_1)
    }
    p_annotation2
  }
