
CorrSamples.Ht.RegionProfile <- function (data, corr.method) {
  ## App fixed terms ##
  dir_funs <- "funs"
  # Preparing objects needed:
  mat <- data
  #########################

  # Preparing corr samples 'cor_samples':
  cor_samples <- cor(mat, method=corr.method) # 'method_cor' is selected by user (pearson, kendall, spearman)
  colnames(cor_samples) <- paste(" -", colnames(cor_samples), "- ", sep="")
  cor_samples_mat <- round(cor_samples, 3) # sample correlations matrix
  #########################

  # Preparing 'p_corr_samples':
  p_corr_samples <- plot_ly(y=rownames(cor_samples),
   x=colnames(cor_samples), z = as.matrix(cor_samples),
    type = "heatmap", text = round(as.matrix(cor_samples),3),
    hoverinfo='x+y+text') %>% colorbar(
    thicknessmode="fraction",
    thickness=0.02,
    ypad=25,
    tickmode="array",
    tickvals=c(-0.96, 0, 0.96),
    ticktext=c(-1, 0, 1),
    ticklen=0
  ) %>% layout(
    font = list(size=10),
    margin = list(l = 120, r = 80, t = 40, b = 40),
    xaxis=list(visible=F)
  )
  #########################
  p_corr_samples
}
