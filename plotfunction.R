suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(ggpubr))


preprocess <- function(ta, meta, methods, phenotype) {
  if(methods != "CIBERSORT") {
    rownames(ta) <- ta[[1]]
    ta <- ta[-1]
    ta <- melt(t(as.matrix(ta)))
  } else {
    ta <- melt(res_ciber)
  }
  
  #combined the phenotype info
  ta <- merge(ta, meta[c("SampleName", phenotype)], by = 1, all = FALSE)
  
  #remove the no expression of immune cell types 
  tmp <- ta %>% group_by(Var2) %>% mutate(mean = mean(value))
  tmp <- tmp[tmp$mean > 0,]
  ta <- tmp[-ncol(tmp)]
  
  ta <- ta[order(ta$Responder, decreasing = FALSE),]
  ta$Var1 <- paste0(ta$Var1, "_",ta[[phenotype]], sep = "")
  ta$Var1 <- factor(ta$Var1, levels = unique(ta$Var1))
  
  return(ta)
}


hmap <- function(ta, meta, methods, phenotype) {
  
  ta <- preprocess(ta, meta, methods, phenotype)
  
  min_value <- 0
  max_value <- round(max(ta$value), digit = 3)
  bar_anno <- c(min_value,max_value)
  
  p <- ggplot(ta, aes_string(ta$Var1, ta$Var2)) +
    geom_point(aes_string(ta$Var1,  ta$Var2, color = ta$value), shape = 15, size = 5) +
    scale_color_gradientn(name = "Score", colours = c("white","firebrick3"))+
    
    theme_light()+theme(axis.text = element_text(colour = 'black', size=10, face = "bold"),
                        axis.title.x = element_blank(),
                        axis.text.x = element_text(angle=90,vjust = 0.5, hjust = 1),
                        axis.title.y = element_blank(),
                        strip.text.x = element_text(size = 10, face = "bold"),
                        legend.title = element_text(size=10, face = "bold"))
  
  return(p)
}

boxfig <- function(ta, meta, methods, phenotype) {
  
  ta <- preprocess(ta, meta, methods, phenotype)
  
  p1 <-  ggplot(ta, aes(x=Var2, y=value, fill = Responder)) +
    geom_boxplot(alpha=0.3, size=0.25,outlier.size= -1, width = 0.5) + theme_bw() +
    stat_compare_means(label = "p.signif", 
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", "")), vjust = 0.5,size = 5) +
    geom_point(position = position_jitterdodge(jitter.width = 0), aes(color = Responder),alpha=0.7) + 
    scale_fill_manual(values=c("#2166AC", "#B2182B")) +
    scale_color_manual(values=c("#2166AC", "#B2182B")) +
    theme(axis.text = element_text(colour = 'black', size=10, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90,vjust = 0.5, hjust = 1),
          #axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title.y = element_text(size = 12, vjust = 2, face = "bold")) + 
    labs(y = paste("Score"))
  
  return(p1)
}


