
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
      palette_score_alpha <- colorRampPalette(c("deepskyblue", "deepskyblue4"), alpha=TRUE)(length(sub_vars))
    } else if (score == "FCS") {
      palette_score_alpha <- colorRampPalette(c("darkorange", "darkorange4"), alpha=TRUE)(length(sub_vars))
    } else if (score == "GCS") {
      palette_score_alpha <- colorRampPalette(c("darkkhaki", "darkolivegreen"), alpha=TRUE)(length(sub_vars))
    }
    
    # p <- ggplot(data, aes(x=reorder(get(var), get(score), FUN=median),
    #                       y=get(score), fill=get(var))) +
    p <- ggplot(data, aes(x=get(var), y=get(score), fill=get(var)))+
      stat_boxplot(geom ='errorbar', width=0.2, color="gray40") +
      # geom_violin(color="gray40", aes(fill=get(var)), trim=FALSE) +
      geom_boxplot(color="gray40", aes(fill=get(var))) +
      scale_fill_manual(values=setNames(palette_score_alpha, totem_order_by_score)) +
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
        
        ans_test <- t.test(values1,values2)
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
      
      p <- p + geom_signif(comparisons = e_list, map_signif_level=FALSE, step_increase=0.06 , annotation=p_vals_here, textsize = 5)
      
    }
    
    
  } else if (ncol(data)==3) {
    if (score == "BCS") {
      color <- "deepskyblue4"
    } else if (score == "FCS") {
      color <- "darkorange4"
    } else if (score == "GCS") {
      color <- "darkolivegreen"
    }
    p <- ggplot(data, aes(x="", y=get(score))) +
      stat_boxplot(geom ='errorbar', width=0.2, color="gray40") +
      # geom_violin(color="gray40", fill=color, trim=FALSE) +      
      geom_boxplot(color="gray40", fill=color) +
      xlab("Global") + ylab(score)
    
  }
  
  
  
  
  p2 <- p + theme_classic()+
    theme(text = element_text(size=14),
          # axis.title.y= element_text(margin=margin(t=0,r=10,b=0,l=0)),
          # axis.title.x= element_text(margin=margin(t=15,r=0,b=0,l=0)),
          axis.text.x = element_text(size=12, angle=45, hjust=1),
          plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="none")

  
  p2
  
}