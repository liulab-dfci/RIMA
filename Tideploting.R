suppressMessages(library(patchwork))
suppressMessages(library(plyr))
suppressMessages(library(optparse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

cmpr_biomk <- function(msi_path, tide_path, meta_path, phenotype) {
  
  #reading the data
  meta <- read.table(meta_path, sep = ",", stringsAsFactors = FALSE, header = TRUE)
  print("Metasheet")
  print(meta)
  msi <- read.table(msi_path, sep="\t", stringsAsFactors = FALSE)
  colnames(msi)<-c("SampleName","MSI")
  print("MSI")
  print(msi)
  
  tide <- read.table(tide_path, sep="\t", stringsAsFactors = FALSE, header = TRUE)
  colnames(tide)[c(1,3)] <- c("SampleName", "Predicted_Response")
  tide <- tide[c(-2,-6)]
  print("TIDE")
  print(tide)
  
  #Combine the MSI score
  tide_table <- merge(tide, msi, by = 1, all = FALSE)
  tide_table <- tide_table[c(1,2,3,7,8,9,10,11,12,4,5,6,13)]
  tide_table$Predicted_Response <- ifelse(tide_table$Predicted_Response == "True", "R","NR")
  
  #whether the data have Responder info
  if ("Responder" %in% colnames(meta)) {
    tide_table <- merge(meta[c("SampleName","Responder")], tide_table, by = 1, all = FALSE)
    colnames(tide_table)[2] <- "Actual_Response"
  }
  
  #Biomarker evaluation plot
  final_table <- merge(tide_table, meta, by = 1, all = FALSE)
  final_table$all <- "z_samples"
  final_table <- rename(final_table, c("MSI" = "MSI Score"))
  
  bioma <- c("MSI Score", "TIDE", "CD8", "CD274")
  new_bioma <- c("Dysfunction", "Exclusion","MDSC","CAF")
  
  plot_list<-lapply(bioma, biomark_process, phenotype)
  plot_list2 <- lapply(new_bioma, biomark_process, phenotype)
  
  combined_plot <- c(plot_list, plot_list2)
  plot <- ggarrange(plotlist = combined_plot, ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom")
  
  return(plot)
}

biomark_process <- function(bioma, phenotype) {
  if ("Timing"%in%phenotype) {
    #preprocess the timing table
    tmp_table <- final_table[c("SampleName", bioma,"PatName", "Responder","Timing", "all")]
    colnames(tmp_table)[2] <- "biomarker"
    tmp_table$Timing <- ifelse(tmp_table$Timing == "pre" | tmp_table$Timing == "Pre", "Apre",
                               ifelse(tmp_table$Timing == "post" | tmp_table$Timing == "Post", "Post", "other"))
    compare_list<-list(unique(as.character(tmp_table[,"Timing"])))
    
    if (length(compare_list[[1]]) > 2) {
      condition <- as.data.frame(compare_list[[1]])
      condition[[1]] <- condition[order(condition[[1]], decreasing = FALSE),]
      
      compare_list <- solution(condition, nrow(condition))
    }
    
    p <- ggplot(tmp_table, aes(x=Timing, y=biomarker),color=Timing) +
      geom_boxplot(aes(fill=factor(Timing)), alpha=0.3, size=0.25,outlier.size= -1, width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = compare_list,aes(label = ..p.format..), tip.length = 0.01,size = 3) +
      geom_point(aes(shape = Responder, color = Responder),alpha=0.7) + geom_line(aes(group = PatName, color = Responder), alpha = 0.5) +
      scale_x_discrete(name ="", labels=c("Apre" = "Pre", "z_samples" = "All Samples")) +
      theme(
        axis.text.x=element_text(size=10, face = "bold", vjust = 0),
        axis.text.y=element_text(size=10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.position="bottom",
        legend.title = element_text(size=10, face = "bold"),
        legend.text = element_text(size=10, face = "bold")
      ) + labs(y = bioma) + guides(fill=FALSE)
  } else {
    tmp_table <- final_table[c("SampleName", bioma,"PatName", phenotype, "all")]
    colnames(tmp_table)[2] <- "biomarker"
    colnames(tmp_table)[4] <- "phenotype"
    compare_list<-list(unique(as.character(tmp_table[, "phenotype"])))
    
    if (length(compare_list[[1]]) > 2) {
      condition <- as.data.frame(compare_list[[1]])
      condition[[1]] <- condition[order(condition[[1]], decreasing = FALSE),]
      
      compare_list <- solution(condition, nrow(condition))
    }
    p <- ggplot(tmp_table, aes(x=phenotype, y=biomarker)) +
      geom_boxplot(aes(fill=factor(phenotype)), alpha=0.3, size=0.25,outlier.size= -1, width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = compare_list,aes(label = ..p.format..), tip.length = 0.01,size = 3) + geom_point(aes(shape = phenotype, color = phenotype),alpha=0.7) +
      geom_boxplot(aes(x = all, y=biomarker, fill=factor(all)), alpha=0.3, size=0.3,outlier.size= -1, width = 0.5)+
      scale_x_discrete(name ="", labels=c("z_samples" = "All Samples")) +
      theme(
        axis.text.x=element_text(size=10, face = "bold", angle = 90, vjust = 0),
        axis.text.y=element_text(size=10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        #legend.position="bottom",
        legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12, face = "bold")
      ) + labs(y = bioma) + guides(fill=FALSE)
  }
}


#function for multiple comparison
solution <- function(condition, len) {
  uni_vector <-  condition[[1]]
  tmp_ta <- NULL
  for (i in c(1:len)) {
    tmp_vector <- uni_vector[i]
    gap <- len - i
    if (gap >= 1) {
      for (n in c(i:(len-1))) {
        tmp_conbination <- data.frame("vector_1" = uni_vector[i], "vector_2" = uni_vector[n+1], "num" = n+1-i)
        tmp_ta <- rbind(tmp_ta, tmp_conbination)
      }
    } else {
      next
    }
  }
  tmp_ta <- tmp_ta[order(tmp_ta$num, decreasing = FALSE),]
  for (i in c(1:nrow(tmp_ta))) {
    tmp <- c(as.character(tmp_ta[i,1]), as.character(tmp_ta[i,2]))
    compare_list[[i]] <- tmp
  }
  return(compare_list)
}


