

survival.plot <- function(data, colors){
  # Survival plot (as percentages [%])
  # - 'data' has 3 columns:  'loss' ; 'normal' ; 'gain'
  #     'rownames(data)' must be regions (or names) to be x-axis
  # - 'colors' is a vector with many colors as columns there are
  # - 'var_surv' is a variable name to define gropus




levels(data[,2]) <- c(1,2) # if Dead/Alive: Alive=1 and Dead=2


# if (colnames(data)[4]=="none") {
#   f1 <- survfit(Surv(time, status==2), data = data)
#   colors <- "darkred"
#   legend_labs <- ""
#   
# } else {
  f1 <- survfit(Surv(time, status==2) ~ variable, data = data)
  legend_labs <- as.character(unique(data[,4]))
  if ( length( which(is.na(legend_labs)) ) > 0 ) {
    legend_labs <- legend_labs[-which(is.na(legend_labs))]
  }
  colors <- colorRampPalette(brewer.pal(5, colors), alpha=TRUE)(length(legend_labs))
  
# }

gg <- ggsurvplot(
  f1, 
  data = data, 
  size = 1,                 # change line size
  palette = colors,# custom color palettes
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.title= "",
  legend.labs = legend_labs,    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme(text = element_text(size=18),
                  legend.position=c(0.9,0.9),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))      # Change ggplot2 theme
) + labs(x="Time")


gg

}