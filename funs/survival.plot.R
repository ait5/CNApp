

survival.plot <- function(data, colors){
  # Survival plot (as percentages [%])
  # - 'data' has 3 columns:  'loss' ; 'normal' ; 'gain'
  #     'rownames(data)' must be regions (or names) to be x-axis
  # - 'colors' is a vector with many colors as columns there are




levels(data[,2]) <- c(1,2) # if Dead/Alive: Alive=1 and Dead=2
legend_labs <- levels(data[,4])

f1 <- survfit(Surv(time, status==2) ~ variable, data = data)


gg <- ggsurvplot(
  f1, 
  data = data, 
  size = 1,                 # change line size
  palette = colors,# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.title= "",
  legend.labs = legend_labs,    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

gg

}