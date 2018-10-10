
Ht.Annotation.RegionProfile <- function(data, mat_vars,
  variables_info_mat, track_vars) {

    # Preparing objects needed:
    mat <- data
    colnames(mat) <- paste(" -", colnames(mat), "- ", sep="")
    sample_names <- colnames(mat)
    seg_regions_names <- rownames(mat)

    var_tracks <- track_vars
    mat_variables <- mat_vars
    mat_annot <- as.matrix(t(mat_variables))

    l_vars_annot_in_plot <-  which(colnames(mat_variables)%in%var_tracks)
    seq_1 <- seq(1,length(l_vars_annot_in_plot), by=1)
    #########################

    plotList <- function(x){
      i <- l_vars_annot_in_plot[x]
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
        if (length(tickvals)>20){
          tickvals <- ceiling(length(tickvals)/2)
          ticktext <- y
        }
      } else if (class_var=="categoric") {
        mat_test <- as.matrix(as.numeric(as.factor(mat_annot[i,])))
        tickvals <- unique(as.numeric(as.character(mat_test)))
        range <- tickvals[length(tickvals)]-tickvals[1]
        tickvals[1] <- tickvals[1] + (range*0.1)
        tickvals[length(tickvals)] <- tickvals[length(tickvals)] - (range*0.1)
        ticktext <- unique(as.character(text))
        if (length(tickvals)>20){
          tickvals <- ceiling(length(tickvals)/2)
          ticktext <- y
        }
      }
      
      if (y == "BCS") {
        palette <- colorRampPalette(c("deepskyblue", "deepskyblue4"), alpha=TRUE)(length(tickvals_0))
      } else if (y == "FCS") {
        palette <- colorRampPalette(c("darkorange", "darkorange4"), alpha=TRUE)(length(tickvals_0))
      } else if (y == "GCS") {
        palette <- colorRampPalette(c("darkkhaki", "darkolivegreen"), alpha=TRUE)(length(tickvals_0))
      } 
      # else if (y == "ID") {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"PiYG"), alpha=TRUE)(length(tickvals))
      # } else if (y == "Lineage") {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"Spectral"), alpha=TRUE)(length(tickvals))
      # } else if ( length(grep("Taylor", y))==1 ) {
      #   colors_groups <- colorRampPalette(brewer.pal(3,"Spectral"), alpha=TRUE)(length(tickvals))
      # }
      

      p1 <- plot_ly(y = y, z = t(mat_test), x = x, type = "heatmap", colors=palette, showscale=T,  text = t(text), hoverinfo='text')%>%
      colorbar(
        thicknessmode="fraction",
        thickness=0.02,

        tickmode="array",
        tickvals= tickvals,
        ticktext= ticktext,
        ticklen=0
      )%>%
      layout(
        font = list(size=10),
        xaxis=list(visible=F),
        yaxis=list(visible=T)
      )

      p1

    }


    jj <- lapply(seq_1, plotList)
    p_annotation <- subplot(jj, nrows=length(jj), shareX=T)
    p_annotation

  }
