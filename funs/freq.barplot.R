

freq.barplot <- function(data, colors){
  # Frequency plot (as percentages [%])
  # - 'data' has 3 columns:  'loss' ; 'normal' ; 'gain'
  #     'rownames(data)' must be regions (or names) to be x-axis
  # - 'colors' is a vector with many colors as columns there are

  mt <- data
  rownames(mt) <- paste(" -", rownames(mt), "- ", sep="")

  x <- as.character(rownames(mt))
  max_term <- max(apply(data, 1, function(x){ sum(as.numeric(as.character(x))) }))
  five_fold_yaxis <- floor(max_term/5)
  tick_vals <- c((max_term - (five_fold_yaxis*5)),
  (max_term - (five_fold_yaxis*4)),
  (max_term - (five_fold_yaxis*3)),
  (max_term - (five_fold_yaxis*2)),
  (max_term - (five_fold_yaxis*1)),
  max_term)
  tick_text <- paste(tick_vals, "  ", sep="")

  p <- plot_ly(mt, x=x, y=mt[,1], name=colnames(mt)[1],type = "bar", marker=list(color=colors[1]), showlegend=F)

  for (i in 2:ncol(mt)){
    p <- p %>% add_trace(y=mt[,i], name=colnames(mt)[i], marker=list(color=colors[i]))
  }
  p <- p %>%   layout(
    yaxis =list(title= "n regions", type="linear",
    range=seq(1:nrow(mt)),
    tickmode="array",
    tickvals=tick_vals,
    ticktext=tick_text
  ),
  xaxis = list(visible=F, categoryarray = x, categoryorder = "array"),
  barmode="stack",

  font = list(size=10),
  margin = list(l = 100, r = 80, t = 40, b = 40)
)
p

}
