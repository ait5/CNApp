############
# Proposta LAIA #
############

### Score.boxplot.by.var.R ###

# score <- "GCS"
# data <- ff_box_score

score.boxplot.by.var <- function (score, data) {
  
  ### SCORE distribution by var:
  # score <- "BCS"
  # print(score)
  
  score_values <- data[,score]
  
  if ( ncol(data)==4) {
    data <- data[order(data[,4]),]
    var <- colnames(data)[4]
    sub_vars <- unique(as.character(data[,var]))
    totem <- reorder(data[,var], data[,score], FUN=median)
    totem_order_by_score <- names(summary(totem))
    
    if (score == "BCS") {
      palette_score_alpha <- colorRampPalette(c("deepskyblue1", "royalblue3"), alpha=TRUE)(length(sub_vars)) ##he canviat els blaus
    } else if (score == "FCS") {
      palette_score_alpha <- colorRampPalette(c("orange1", "orangered2"), alpha=TRUE)(length(sub_vars)) ##he canviat els taronges
    } else if (score == "GCS") {
      palette_score_alpha <- colorRampPalette(c("limegreen", "springgreen4"), alpha=TRUE)(length(sub_vars)) ##he canviat els verds
    }
    
    # p <- ggplot(data, aes(x=reorder(get(var), get(score), FUN=median),
    #                       y=get(score), fill=get(var))) +
    p <- ggplot(data, aes(x=get(var), y=get(score), fill=get(var)))+
      stat_boxplot(geom ='errorbar', width=0.2, color="black") +    ##canvi gray40 to black
      # geom_violin(color="gray40", aes(fill=get(var)), trim=FALSE) +
      geom_boxplot(color="black", size = 1, width = 0.6, aes(fill=get(var))) +   ##canvi gray40 to black and add size / width parameters
      scale_fill_manual(var, values=setNames(palette_score_alpha, totem_order_by_score)) +
      xlab(var) + ylab(score)
    
    real_groups <- sub_vars
    
    if (length(real_groups)<=6) {
      
      e_list <- list()
      for (i in 1:length(real_groups)){
        for (j in (i:length(real_groups))[-(i:length(real_groups)==i)]){
          contr <- c(as.character(real_groups[i]),as.character(real_groups[j]))
          e_list <- c(e_list, list(contr))
        }
      }
      
      p_vals_here <- as.vector(do.call(cbind,lapply(e_list, function(x){
        name1 <- x[1]
        name2 <- x[2]
        values1 <- data[which(data[,var]==name1),score]
        values2 <- data[which(data[,var]==name2),score]   
        
        ans_test <- wilcox.test(values1,values2) # Wilcox-test
        p_value <- round(ans_test$p.value,3)
        p_value
      })))
      
      
      for (i in 1:length(p_vals_here)){
        val <- as.numeric(as.character(p_vals_here[i]))
        if (val <= 0.001){
          val2 <- paste("***", val, sep=" ")
        } else if (val <= 0.01){
          val2 <- paste("**", val, sep=" ")
        } else if (val <= 0.05){
          val2 <- paste("*", val, sep=" ")
        } else if (val > 0.05){
          val2 <- paste("(ns)", val, sep=" ")
        } else {
          val2 <- val
        }
        p_vals_here[i] <- val2
      }
      
      p <- p + geom_signif(comparisons = e_list, map_signif_level=FALSE, step_increase=0.06 , annotation=p_vals_here, textsize = 5.5)
      
    }
    
    
  } else if (ncol(data)==3) {
    if (score == "BCS") {
      color <- "deepskyblue1"          ##canvi de color (com el que he posat a dalt)
      color_line <- "royalblue3"
    } else if (score == "FCS") {
      color <- "orange1"          ##canvi de color (com el que he posat a dalt)
      color_line <- "orangered2"
    } else if (score == "GCS") {
      color <- "limegreen"        ##canvi de color (com el que he posat a dalt)   
      color_line <- "springgreen4"
    }
    # p <- ggplot(data, aes(x="", y=get(score))) +
    #   stat_boxplot(geom ='errorbar', width=0.2, color="black") + ##canvi gray40 a black
    #   # geom_violin(color="gray40", fill=color, trim=FALSE) +      
    #   geom_boxplot(color="gray40", fill=color,  size = 1, width = 0.6) +   ##afegit size and width
    #   xlab("Global") + ylab(score)
    v_median <- median(data[,score])
    
    p <- ggplot(data, aes(x=get(score))) +
      geom_density(color="black", size = 1, width = 0.6, fill=color)+
      geom_vline(xintercept = v_median, color=color_line, linetype="dashed", size=0.6)+
      theme(legend.position="none")+
      labs(x=paste(score, "values"))+
      theme(text = element_text(size=40),
            legend.position=c(0.8,0.8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
    
    
  }

  
  #new version
  p2 <- p + theme_minimal()                  ##canvi de theme
  theme(axis.text.x = element_text(colour="black", size = 14, face='bold'), ##canvi linies / lletra
        axis.text.y = element_text(colour="black", size = 14), 
        axis.title.y = element_text(size = 14, face = 'bold.italic'),
        legend.position="none")
  
  
  
  p2
  
}
