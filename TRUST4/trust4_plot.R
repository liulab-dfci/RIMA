suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(reshape))
suppressMessages(library(ggpubr))
suppressMessages(library(GGally))
suppressMessages(library(network))
suppressMessages(library(ggcorrplot))
suppressMessages(library(sna))

CompareGroups <- function(dat,F1,F2,cols,y_pos,metric){
  par(oma=c(2,2,2,2))
  gp <- ggplot(dat, aes_string(x = F1, y = F2, fill = F1))+
    #geom_violin(trim = TRUE,alpha = 0.3,scale = "width",lwd=0.2,width=0.3)+
    geom_boxplot(lwd=0.3,size=0.3,outlier.size=-1,width=0.3,alpha=0.3,position = position_dodge2(preserve = "single"))+
    theme_bw() +
    ylab(metric)+
    scale_fill_manual(values = brewer.pal(n = 8, name = "Set2")[1:cols])+
    stat_compare_means( aes(label = ..p.format..),
                        label.x = 1.4, label.y = y_pos,
                        tip.length = 0.01,size = 2)+
    scale_x_discrete(name ="", labels=c("Apre" = "Pre")) +
    
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=10,r=5,t=4,b=4),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70,size=8,face = "bold",hjust=1),
          axis.text.y=element_text(size=8,hjust=1,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8,face="bold"),
          legend.position = "bottom",
          legend.title = element_text(size=8, face = "bold"),
          legend.text = element_text(size=8, face = "bold")) + guides(fill=FALSE)
  theme(plot.margin = unit(c(0,0.05,0,0.05), "cm"))
  return(gp)
}


Trust4_plot <- function(phenotype, metasheet) {
  
  clinic.col <- phenotype
  
  cols <- length(unique(meta[,clinic.col]))
  
  
  ##############-------------Metric: fraction of BCR reads------------------------###############
  ################################################################################################
  cdr3.lib <- merge(bcr.lib.reads,meta,by.x='sample',by.y='SampleName')
  y_pos <- max(cdr3.lib$Infil)*1.2
  #y_pos<- 0.0001
  gr <- CompareGroups(cdr3.lib,clinic.col,"Infil",cols,y_pos,"Fraction of BCR reads")
  if (nrow(bcr.lib.reads) == 0) {
    gr <- NULL
  }
  
  if(clinic.col=="Timing") {
    
    gr <- gr + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) + 
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
  }
  
  ##############-------------Metric: unique BCR cdr3s size ------------------------###############
  ################################################################################################
  tmp <- aggregate(CDR3aa ~ sample, cdr3.bcr.heavy, function(x) length(unique(x))) 
  cdr3.size <- merge(tmp,meta,by.x='sample',by.y='SampleName')
  #y_pos=100
  y_pos <- max(cdr3.size$CDR3aa)*1.2
  gc <- CompareGroups(cdr3.size,clinic.col,"CDR3aa",cols,y_pos,"Unique CDR3")
  
  if(clinic.col=="Timing") {
    gc <- gc + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
  }
  
  ##############-------------Metric: somatic hypermutation rate ------------------------###########
  ################################################################################################
  tmp <- na.omit(cbind.data.frame(sample = names(SHMRatio), SHM.rate = unlist(SHMRatio)))
  tmp['name'] <-   str_replace(tmp$sample, "_TRUST4_BCR_heavy_cluster.Rdata", "")
  
  all.SHM <- merge(tmp,meta,by.x='name',by.y='SampleName')
  
  y_pos <- max(all.SHM$SHM.rate)*1.0
  gs <- CompareGroups(all.SHM,clinic.col,"SHM.rate",cols,y_pos,"SHM rate")
  
  if(clinic.col=="Timing") {
    gs <- gs + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
  }
  
  ##############-------------Metric: BCR clonality ------------------------###########################
  ################################################################################################
  tmp <- data.frame(do.call("rbind",bcr_clonality)) %>% 
    mutate(clonality = as.numeric(as.character(clonality)))
  
  clonality.all <- merge(tmp,meta,by.x='sample',by.y='SampleName')
  
  if ("NaN"%in%clonality.all$clonality) {
    clonality.all <- subset(clonality.all, clonality.all[[2]] != "NaN")
    
  }
  
  
  y_pos <- max(clonality.all$clonality)*1.0
  gn <- CompareGroups(clonality.all,clinic.col,"clonality",cols,y_pos,"Clonality")
  #gn
  
  if(clinic.col=="Timing") {
    gn <- gn + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size = 0.3)
  }
  
  bcr_plot <- ggarrange(gr, gc, gs, gn, ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
  
  
  ##############-------------Metric: fraction of TCR reads ------------------------###########################
  ################################################################################################
  cdr3.lib <- tcr.lib.reads %>% mutate(clinic = meta[sample,clinic.col])
  y_pos <- max(cdr3.lib$Infil)+0.0000001
  
  if (clinic.col == "Timing") {
    cdr3.lib <- merge(cdr3.lib, meta[c("SampleName", "Responder", "PatName")], by = 1, all = FALSE)
  }
  gr <- CompareGroups(cdr3.lib,"clinic","Infil",cols,y_pos,"Fraction of TCR reads")
  if (nrow(tcr.lib.reads) == 0) {
    gr <- NULL
  }
  
  if(clinic.col=="Timing") {
    gr <- gr + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
  }
  
  ##############-------------Metric: TCR unique cdr3 size ------------------------###########################
  ################################################################################################
  tmp <- aggregate(CDR3aa ~ sample, cdr3.tcr, function(x) length(unique(x))) 
  cdr3.size <- merge(tmp,meta,by.x='sample',by.y='SampleName') 
  
  y_pos <- max(cdr3.size$CDR3aa)-0.4
  gc <- CompareGroups(cdr3.size,clinic.col,"CDR3aa",cols,y_pos,"Unique CDR3")
  
  if(clinic.col=="Timing") {
    gc <- gc + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
  }
  
  ##############-------------Metric: TCR CPK ------------------------###########################
  ################################################################################################
  tmp <- aggregate(CDR3aa ~ sample+lib.size, cdr3.tcr, function(x) length(unique(x))) %>%
    mutate(CPK = signif(CDR3aa/(lib.size/1000),4))
  
  cpk <-  merge(tmp,meta,by.x='sample',by.y='SampleName') 
  y_pos <- max(cpk$CPK)-0.4
  
  gk <- CompareGroups(cpk,clinic.col,"CPK",cols,y_pos,"Clonetypes per kilo reads")
  
  if(clinic.col=="Timing") {
    gk <- gk + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
  }
  
  ##############-------------Metric: TCR clonality ------------------------###########################
  ################################################################################################
  tmp <- data.frame(do.call("rbind",tcr_clonality)) %>% 
    mutate(clonality = as.numeric(as.character(clonality)))
  
  clonality.all <- merge(tmp,meta,by.x='sample',by.y='SampleName')
  if ("NaN"%in%clonality.all$clonality) {
    clonality.all <- subset(clonality.all, clonality.all[[2]] != "NaN")
  }
  y_pos <- max(clonality.all$clonality)*1.1
  
  gn <- CompareGroups(clonality.all,clinic.col,"clonality",cols,y_pos,"Clonality")
  
  if(clinic.col=="Timing") {
    gn <- gn + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
      geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
  }
  
  tcr_plot <- ggarrange(gr, gc, gk, gn,ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
  
  p <- ggarrange(bcr_plot, tcr_plot, nrow = 2)
  
  return(p)
}




