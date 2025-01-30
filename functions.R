

################################################################################
#                          ALPHA DIVERSITY                                     #
################################################################################

# This function calculate alpha diversity metrics for a phyloseq object
# 
# phy_objetc <- Phyloseq object created with qiime2 data
# group <- metadata column you want to use to separate data.
# alfa_metrics <- vector with the statistical metrics to calculate alfa diversity
# my comparision <- a list of pairs of condition you want to compare
# stat_metrics <- statistica metric do you want to use to see if there is a significant difference between conditions
# asterisk <- if yes show astherisks in accordance to significance level, if no write the pvalue

alpha_diversity_plot <- function(phy_object, group, alfa_metrics, my_comparisons, stat_metric, asterisk, title){
  library(ggpubr)
  #alpha diversity plot
  p <- plot_richness(phy_object, x=group, measures=alfa_metrics) 
  p <- p + geom_boxplot(alpha=9)
  p <- p + theme_light()
  
  if (asterisk == TRUE){
    p <- p + stat_compare_means(method = stat_metrics, comparisons = my_comparisons, label = "p.signif",  p.adjust.method = "BH")
  }else{
    p <- p + stat_compare_means(method = stat_metrics, p.adjust.method = "BH")  
    
  }
  p <- p + ggtitle(label = title) + theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title=element_text(size=14,face="bold"),
    axis.text=element_text(size=12),
    strip.text = element_text(size = 14),
    axis.title.x=element_blank())
  
  return(p)
}
  

################################################################################
#                           ABUNDANCE                                          #
################################################################################

plot_abundance_by_taxa <- function(df, taxa, metric, group){
  plt <- ggplot(df ,aes(x=get(taxa), y=get(metric), group=get(group), fill=get(group))) +
    geom_bar(colour="black", stat="identity",position="dodge")  +
    theme_classic() +
    theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1, face = "italic"), 
          legend.title = element_blank(),
          axis.title=element_text(size=20,face="bold"),
          axis.text.y=element_text(size=18),
          legend.text=element_text(size=20),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_fill_grey(start = 0, end = 1) +
    geom_errorbar(aes(ymin=get(metric), ymax=get(metric)+sd),
                  width=.5,                    # Width of the error bars
                  position=position_dodge(.9)) 
  
  return(plt)
  
}

plot_abundance_heatmap<- function(df, metric, taxa_level, group){
  ggplot(df, aes(get(group), get(taxa_level), fill= get(metric))) + 
    geom_tile() +
    scale_fill_distiller(palette = "Spectral", name = metric) +
    ylab(taxa_level)+
    theme(panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x  = element_text(angle=90, vjust=0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))
}

################################################################################
#                          DIFFERENTIAL ABUNDANCE DESEQ2                       #
################################################################################


daa_bar_plot <- function(df, contrast, title){
  theme_set(theme_bw())
  
    #Transform data
  # Phylum order
  x = tapply(df$log2FoldChange, df$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  df$Phylum = factor(as.character(df$Phylum), levels=names(x))
  # Genus order
  x = tapply(df$log2FoldChange, df$Genus, function(x) max(x))
  x = sort(x, TRUE)
  df$Genus = factor(as.character(df$Genus), levels=names(x))
  
    #Generate plot
  ggplot(df) +
    geom_col(aes(x = log2FoldChange, y = Genus, fill = Phylum)) + 
    geom_vline(xintercept = 0.0, color = "Black", size = 0.7)  +
    ggtitle(title) +
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.title=element_text(size=14,face="bold"),
      axis.text=element_text(size=12),
      legend.text=element_text(size=12)) +
    theme_minimal()
}

################################################################################
#                          DIFFERENTIAL ABUNDANCE                              #
################################################################################

lefse_lda_plot <- function(df, title){
 ggplot(df) +
 geom_bar(aes(x = ef_lda, y = feature, fill = enrich_group), stat="identity",position="dodge", width = 0.5) +
 xlab("LDA score (log10)") + ylab(" ") +
 ggtitle(title) + 
 theme(legend.title=element_blank())
}


################################################################################
#                          BETA DIVERSITY                                      #
################################################################################

betadiversity_analysis <- function(pseq, method, distance, weighted = FALSE){
  if(tolower(distance) == "unifrac"){
    if(weighted == TRUE){
      phyloseq::ordinate(pseq, method = method, distance = distance, weighted=T)
    }else{
      phyloseq::ordinate(pseq, method = method, distance = distance, weighted=F)
    }
  }else{
    if(weighted == FALSE){
      phyloseq::ordinate(pseq, method = method, distance = distance)
    }else{
      print("These metrics do not use weight")
    }
  }
}

plot_distance <- function(pseq, bd, color, shape=NULL, elipse = TRUE){
  if(is.null(shape) == FALSE){
    plt <- plot_ordination(pseq, bd, color= color, shape=shape) + geom_point(size=3)
  }else{
    plt <- plot_ordination(pseq, bd, color= color) + geom_point(size=3)
  }
  
  plt <- plot_ordination(pseq, bd, color= color, shape=shape) + geom_point(size=3)
  if (elipse ==TRUE){
    plt <- plt + stat_ellipse() 
  }
  
  plt <- plt + theme_classic() +theme(axis.title=element_text(size=14,face="bold"),
                     axis.text=element_text(size=12),
                     legend.text=element_text(size=12),
                     legend.title=element_text(size=14))
  return(plt)
}
