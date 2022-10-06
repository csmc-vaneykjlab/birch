library(shiny)
library(shinyjs)
library(shinyalert)
library(shinyvalidate)
library(shinythemes)
library(shinycssloaders)
library(rsconnect)
library(DT)
library(tidyr)
library(textreadr)
#library(rjson)
library(optparse)
library(dplyr)
library(tibble)
library(ggplot2)
library(assertthat)
library(readr)
library(purrr)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(DEP)
library(proBatch)
library(gtable)
library(grid)
library(gridExtra)
#library(sqldf)
library(pcaMethods)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(tools)
library(missRanger)
#library(iheatmapr)
library(plotly)

############## function to help with validation
`%then%` <- function(a, b) {
  if (is.null(a)) b else a
}

########## Functions from bc script
#Replace hyphens in colnames with underscore
colClean <- function(x){ 
  colnames(x) <- gsub("-", "_", colnames(x))
  colnames(x) <- gsub("\\.", "_", colnames(x))
  return(x)  
}

# Replace hyphens in Level3 of annotation with underscore
colClean2 <- function(anno_dataframe, sample_names){ 
  levels(anno_dataframe[[sample_names]]) <- gsub("-", "_", levels(anno_dataframe[[sample_names]]))
  levels(anno_dataframe[[sample_names]]) <- gsub("\\.", "_", levels(anno_dataframe[[sample_names]]))
  #levels(x$Level3) <- gsub("-", "_", levels(x$Level3))
  #levels(x$Level3) <- gsub("\\.", "_", levels(x$Level3))
  return(anno_dataframe)
}

# rename_level3, check_headers, filter_values, protein_concat, order_samples
#rename Level3 in annotation.txt to FullRunName
#NOTE: Needs to be called before filter_values function
rename_level3 <- function(anno_dataframe, sample_names){
  names(anno_dataframe)[names(anno_dataframe) == sample_names] <- "FullRunName"
  return(anno_dataframe)
}

# Create color_list 
create_color_list <- function(annotation, cols_of_int, exp_grp, sample_names) {
  check_cols_of_interest <- unlist(strsplit(cols_of_int, " "))
  
  suppressMessages({
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors,technical_factors)
    plot_for_pca <- c(biological_factors,technical_factors)
    
    color_list <- sample_annotation_to_colors(
      annotation,
      factor_columns = c(biological_factors,technical_factors),
      numeric_columns = c('order')
    )
  })
  return(color_list)
}

# Function to create bar graphs for initial analysis
plot_bar <- function(annotation_data, col_of_int, color_list) {
  # necessary packages / libraries
  require(ggplot2)
  require(RColorBrewer)
  require(ggpubr)
  
  table1_stacked <- annotation_data %>% 
    group_by("Batch_name"=annotation_data[[col_of_int]]) %>%
    dplyr::tally() %>%
    dplyr::mutate(percent=n/sum(n)) %>%
    mutate(Batch_type = toString(col_of_int))
  
  p <- ggplot(data = table1_stacked, aes(x = Batch_type, y=n, fill = Batch_name)) +
    geom_bar(stat="identity", position = "fill", width=0.4) +
    geom_label(aes(label = paste0(sprintf("%1.1f", percent*100),"%, n=", n), group = Batch_name), position = position_fill(vjust = 0.5), fill="white", size=2) +
    scale_fill_manual(values = color_list[[col_of_int]]) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_minimal()
  
  p <- p + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          #legend.position="bottom",
          legend.position = c(0.8, 0.2),
          legend.title=element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    ylab("Percent") +
    labs(fill = toString(col_of_int)) 
  
  return(p)
}

# function to process NAs in the data
call_missing_vals <- function(sample_matrix, annotation_data, group_by, protein) {
  results_matrix <- sample_matrix
  results_matrix[is.na(results_matrix)] = 0
  
  overlaps <- pivot_longer(results_matrix, cols = c(everything(), -protein), names_to = "FullRunName", values_to = "Intensity") %>% 
    left_join(annotation_data, by="FullRunName") %>% 
    dplyr::rename("SampleName" = "FullRunName" , "Group" := !!group_by ) %>% 
    subset(select=c(protein,"Intensity","SampleName","Group"))
  return(overlaps)
}

# Balloon plot showing a grid of exp grp X batch col  
plot_balloon <- function(annotation_data, batch_col, exp_grp) {
  anno_for_plot <- annotation_data %>%  dplyr::count(across(exp_grp), across(batch_col))
  
  p <- ggplot(anno_for_plot, aes(x=anno_for_plot[[exp_grp]], y=anno_for_plot[[batch_col]])) +
    geom_point(aes(size = n), shape = 21, colour = "black", fill = "cornsilk") +
    scale_size_area(max_size = 15, guide = "none") +
    geom_text(aes(
      y = as.numeric(as.factor(anno_for_plot[[batch_col]])) - sqrt(n)/34, label = n),
      vjust = 2.0,
      colour = "black",
      size = 3
    ) +
    ggtitle(paste0("Experimental group X ", batch_col)) +
    xlab(toString(exp_grp)) +
    ylab(toString(batch_col)) +
    theme(plot.title = element_text(size = 15, hjust = 0.5)) +
    theme(plot.margin = margin(1,1,1.5,1, "cm"))
  return(p)
}

# function to make pareto plot
plot_pareto <- function(overlaps, protein) {
  # Get percent missingness per fragment
  missval <- overlaps %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::group_by(across(protein)) %>%
    dplyr::summarize(count = n()-sum(Intensity), percent = (1-(sum(Intensity)/n()))*100)
  
  # Bin missingness and get frequency and cumulative frequency
  percent_miss <- missval$percent
  bins <- seq(0,100,by=10)
  scores <- cut(percent_miss,bins,include.lowest=TRUE)
  freq_table <- transform(table(scores))
  cum_table <- transform(freq_table,Cum_Freq=cumsum(Freq))
  
  # Make pareto plot
  scaleRight <- tail(cum_table$Cum_Freq, n=1)/head(cum_table$Freq, n=1)
  
  ggplot(cum_table, aes(x=cum_table$score)) +
    geom_bar(aes(y=cum_table$Freq), fill='blue', stat="identity") +
    geom_point(aes(y=cum_table$Cum_Freq), color = rgb(0, 1, 0), pch=16, size=1) +
    geom_path(aes(y=cum_table$Cum_Freq, group=1), colour="red", lty=3, size=1) +
    scale_y_continuous(name = "Fragment missingness count", sec.axis = sec_axis(~./ max(cum_table$Cum_Freq), labels = scales::percent, name="Fragment missingness percentage")) +
    theme(text = element_text(size=14, face="bold"), axis.text.x = element_text(angle=90, vjust=0.5, hjust = 0.5, size = 10)) +
    labs(title = "Missingness Distribution", x = 'Missingness distribution range in percentage', y ='Fragment missingness count', size = 12, face = "bold", vjust=0.5, hjust = 0.5)
}

# function to make heatmap showing missing data
plot_missval <- function(overlaps, FONTSIZE_sample_axis, color_list, batch_col, protein) {
  require(tidyverse)
  require(ComplexHeatmap)
  require(RColorBrewer) 
  
  colors <- color_list[[batch_col]]
  colors_df <- as.data.frame(colors)
  colors_df <- colors_df %>%
    rownames_to_column()
  
  n_groups <- length(unique(overlaps$Group))
  GROUP_COLORS <- as.vector(colors_df$colors)
  names(GROUP_COLORS) <- unique(colors_df$rowname)
  
  #GROUP_COLORS <- colorRampPalette(brewer.pal(min(8, n_groups), "Dark2"))(n_groups)
  #names(GROUP_COLORS) <- unique(overlaps$Group)
  
  annot_group <- overlaps %>% 
    dplyr::distinct(SampleName, Group) %>% 
    dplyr::arrange(Group) %>% 
    dplyr::pull(Group)
  
  missval <- overlaps %>%
    dplyr::group_by(across(all_of(protein))) %>%
    mutate(Intensity = ifelse(Intensity > 0, TRUE, FALSE)) %>% 
    dplyr::filter(any(Intensity == F)) %>%
    dplyr::arrange(Group) %>%
    dplyr::select(-Group) %>%
    dplyr::mutate(Intensity = if_else(Intensity, 1, 0)) %>%
    pivot_wider(names_from = SampleName, values_from = Intensity) %>%
    column_to_rownames(var=as.character(protein)) %>%
    as.matrix()
  
  # missval[,order(colnames(missval))]
  Heatmap(
    missval,
    col = c("#FFFFFF", "#000000"),
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_title = "Fragments missing\nfrom at least one sample",
    column_title = paste0("Missing Fragments Pattern: ", batch_col),
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_names_centered = FALSE,
    column_names_gp = gpar(fontsize = FONTSIZE_sample_axis),
    name = "Legend",
    heatmap_legend_param = list(labels = c("Missing", "Observed")),
    cluster_rows = F,
    cluster_columns = F,
    bottom_annotation = HeatmapAnnotation(
      Group = annot_group,
      col = list(Group = GROUP_COLORS)
    )
  )
}

# functions to make fragment filtering plots
plot_missgroup <- function(overlaps, expgrp_threshold, protein) {
  require(tidyverse)
  require(ggpubr)
  
  # Get total counts for each exp group
  annot_group <- overlaps %>% 
    dplyr::distinct(SampleName, Group) %>%
    group_by(Group) %>% 
    mutate(count = n()) %>%
    dplyr::select(-SampleName) %>%
    dplyr::distinct(Group, count) 
  
  
  #rename protein column to "Protein"
  names(overlaps)[names(overlaps) == protein] <- "Protein"
  
  # Get counts of each exp group per fragment
  missval <- overlaps %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::arrange(Group) %>%
    dplyr::select(-SampleName) %>%
    group_by(Protein, Group) %>% 
    summarise(Intensity = sum(Intensity)) %>% 
    left_join(annot_group, by="Group") %>% 
    mutate(status = 1-(Intensity/count)) %>%
    dplyr::select(Protein,status) %>%
    group_by(Protein) %>%
    mutate(status = ifelse(all(status > expgrp_threshold), "Threshold failed", "Threshold passed")) %>%
    dplyr::distinct(Protein,status) %>%
    dplyr::select(-Protein) %>%
    group_by(status) %>% 
    mutate(count = n()/nrow(.)) %>% 
    dplyr::distinct(status,count)
  
  ggplot(missval, aes(x = "Group", y = count, fill = status, label = status)) +
    geom_bar(stat = "identity") +
    labs(x = NULL, y = "% Fragments") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    theme(text = element_text(size=14, face="bold"), axis.text.x = element_text(angle=90, vjust=0.5, hjust = 0.5, size = 14)) +
    ggtitle("Group Missingness") 
}

plot_missbatch <- function(overlaps, batch_threshold, protein) {
  require(tidyverse)
  require(ggpubr)
  
  # Get total counts for each exp group
  annot_group <- overlaps %>% 
    dplyr::distinct(SampleName, Group) %>%
    group_by(Group) %>% 
    mutate(count = n()) %>%
    dplyr::select(-SampleName) %>%
    dplyr::distinct(Group, count) 
  
  #rename protein column to "Protein"
  names(overlaps)[names(overlaps) == protein] <- "Protein"
  
  # Get counts of each exp group per fragment
  missval <- overlaps %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::arrange(Group) %>%
    dplyr::select(-SampleName) %>%
    group_by(Protein, Group) %>% 
    summarise(Intensity = sum(Intensity)) %>% 
    left_join(annot_group, by="Group") %>% 
    mutate(status = 1-(Intensity/count)) %>%
    dplyr::select(Group,status) %>%
    #group_by(Protein) %>%
    mutate(status = ifelse((status > batch_threshold), "Threshold failed", "Threshold passed")) %>%
    #dplyr::distinct(Protein,status) %>%
    dplyr::select(-Protein) %>%
    group_by(Group,status) %>% 
    mutate(count = n()) %>% 
    dplyr::distinct(status,count) %>% 
    group_by(Group) %>%
    mutate(percent = count/sum(count)*100) %>%
    ungroup() %>%
    dplyr::select(-count)
  
  # Get cumulative counts 
  cumulative <- overlaps %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::arrange(Group) %>%
    dplyr::select(-SampleName) %>%
    group_by(Protein, Group) %>% 
    summarise(Intensity = sum(Intensity)) %>% 
    left_join(annot_group, by="Group") %>% 
    mutate(status = 1-(Intensity/count)) %>%
    dplyr::select(Protein,status) %>%
    group_by(Protein) %>%
    mutate(status = ifelse(any(status > batch_threshold), "Threshold failed", "Threshold passed")) %>%
    dplyr::distinct(Protein,status) %>%
    dplyr::select(-Protein) %>%
    group_by(status) %>% 
    mutate(percent = n()/nrow(.)*100) %>% 
    dplyr::distinct(status,percent)
  
  cumulative <- data.frame(append(cumulative, c(Group='Total'), after=1))
  
  missval = bind_rows(missval, cumulative)
  
  ggplot(missval, aes(x = Group, y = percent, fill = status, label = percent)) +
    geom_bar(stat = "identity") +
    labs(x = "Batch", y = "% Fragments") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    theme(text = element_text(size=12, face="bold"), axis.text.x = element_text(angle=90, vjust=0.5, hjust = 0.5, size = 10)) +
    ggtitle("Batch Missingness")
  
}

plot_per_sample <- function(overlaps, FONTSIZE_sample_axis, SAMPLE_NA_CUTOFF, color_list, batch_col) {
  require(tidyverse)
  require(ggpubr)
  require(RColorBrewer)
  
  colors <- color_list[[batch_col]]
  colors_df <- as.data.frame(colors)
  colors_df <- colors_df %>%
    rownames_to_column()
  
  n_groups <- length(unique(overlaps$Group))
  GROUP_COLORS <- as.vector(colors_df$colors)
  names(GROUP_COLORS) <- unique(colors_df$rowname)
  
  #n_groups <- length(unique(overlaps$Group))
  #GROUP_COLORS <- colorRampPalette(brewer.pal(min(8, n_groups), "Dark2"))(n_groups)
  #names(GROUP_COLORS) <- unique(overlaps$Group)
  
  stat2 <- overlaps %>% group_by(SampleName, Group) %>%
    mutate(Intensity = ifelse(Intensity > 0, TRUE, FALSE)) %>% 
    dplyr::summarize(n = n(), sum = (n()-sum(Intensity))/n()*100, .groups = "drop")
  
  ggplot(stat2, aes(x = SampleName, y = sum, fill = Group, order = Group)) +
    geom_col() +
    geom_hline(aes(yintercept = SAMPLE_NA_CUTOFF*100), colour="#990000", linetype="dashed") + 
    scale_fill_manual(values = GROUP_COLORS) +
    labs(title = "Sample Missingness", x = "Sample Name",
         y = "% Missingness") +
    theme_classic2() +
    theme(plot.title = element_text(size=14,face="bold"),
          axis.text.x = element_text(size=FONTSIZE_sample_axis,face="bold", vjust=0.5,angle=90)) +
    coord_cartesian(ylim = c(0, 100))
}

# function for minimal filtering
minimal_filtering <- function(sample_matrix, annotation_data, group_by, protein) {
  results_matrix <- sample_matrix
  results_matrix[is.na(results_matrix)] = 0
  
  overlaps <- pivot_longer(results_matrix, cols = c(everything(), -protein), names_to = "FullRunName", values_to = "Intensity") %>% 
    left_join(annotation_data, by="FullRunName") %>% 
    dplyr::rename("SampleName" = "FullRunName" , "Group" := !!group_by ) %>% 
    subset(select=c(protein,"Intensity","SampleName","Group"))
  
  #rename protein column to "Protein"
  names(overlaps)[names(overlaps) == protein] <- "Protein"
  
  missval <- overlaps %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::arrange(Group) %>%
    dplyr::select(-SampleName) %>%
    group_by(Protein, Group) %>% 
    summarise(Intensity = sum(Intensity)) %>% 
    mutate(check = ifelse(Intensity < 2, 0, 1)) %>% 
    group_by(Protein) %>% 
    filter(any(check == 0)) %>%
    summarise(Intensity = sum(Intensity))
  
  # missval now contains only those proteins/fragments that need to be removed because they have batches with < 2 intensities in them
  rownames(missval) <- missval$Protein
  frags_to_remove = row.names(missval)
  # removing from matrix if missval contains the fragment
  sample_matrix = sample_matrix[!(sample_matrix$Protein %in%  frags_to_remove),]
  
  return(sample_matrix)
}

# function to filter fragments using thresholds
filter_samples <- function(df, annotation, sample_missingness, expgroup_missingness, batch_type, batch_missingness, exp_grp, protein) {
  
  #sprintf("Initial samples = %i, fragments = %i", nrow(annotation), nrow(df))
  
  print("Initial samples:")
  print(nrow(annotation))
  print("Initial fragemnts:")
  print(nrow(df))
  print("\n")
  
  ####### A #######
  samples <- list(annotation$FullRunName) 
  samples = unlist(samples, use.names=FALSE)
  
  sam_remove = df[ lapply( df, function(x) sum(is.na(x)) / length(x) ) >= sample_missingness ]
  samples_to_remove = names(sam_remove)
  
  # drop columns based on names from the original dataframe 
  fil_df = dplyr::select(df, -samples_to_remove)
  
  # remove rows from copy of annotation as well
  annotation = annotation[!(annotation$FullRunName %in% samples_to_remove),]
  row.names(annotation) <- NULL 
  
  ##### LOGGING ######
  
  #sprintf("Drops from a) = %i samples, Remaining = %i samples, %i fragments", length(samples_to_remove), nrow(annotation), nrow(fil_df))
  
  print("Drops from a)")
  print("Samples dropped: ")
  print(length(samples_to_remove))
  print("Remaining Samples:")
  print(nrow(annotation))
  print("Remaining Fragments:")
  print(nrow(fil_df))
  print("\n")
  
  ####### B #######
  
  overlaps2 <- call_missing_vals(fil_df, annotation, exp_grp, protein)
  
  #rename protein column to "Protein"
  names(overlaps2)[names(overlaps2) == protein] <- "Protein"
  
  annot_group <- overlaps2 %>% 
    dplyr::distinct(SampleName, Group) %>%
    group_by(Group) %>% 
    mutate(count = n()) %>%
    dplyr::select(-SampleName) %>%
    dplyr::distinct(Group, count) 
  
  # Get counts of each exp group per fragment
  missval <- overlaps2 %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::arrange(Group) %>%
    dplyr::select(-SampleName) %>%
    group_by(Protein, Group) %>% 
    summarise(Intensity = sum(Intensity)) %>% 
    left_join(annot_group, by="Group") %>% 
    mutate(status = 1-(Intensity/count)) %>%
    dplyr::select(Protein,status) %>%
    mutate(check = ifelse(status <= expgroup_missingness, 1, 0)) %>%
    group_by(Protein) %>% 
    summarise(check = sum(check))
  
  dropped = length(which(missval$check==0))
  
  rownames(missval) <- missval$Protein
  frags_to_keep_c = row.names(missval)[which(missval$check!=0)]
  fil_df = fil_df[(fil_df[[protein]] %in% frags_to_keep_c),]
  
  ##### LOGGING ######
  
  #sprintf("Drops from b) = %i fragments, Remaining = %i samples, %i fragments", dropped, nrow(annotation), nrow(fil_df))
  
  print("Drops from b)")
  print("Fragments dropped: ")
  print(dropped)
  print("Remaining Samples:")
  print(nrow(annotation))
  print("Remaining Fragments:")
  print(nrow(fil_df))
  print("\n")
  
  ####### C #######
  overlaps <- call_missing_vals(fil_df, annotation, batch_type, protein)
  
  #rename protein column to "Protein"
  names(overlaps)[names(overlaps) == protein] <- "Protein"
  
  annot_group <- overlaps %>% 
    dplyr::distinct(SampleName, Group) %>%
    group_by(Group) %>% 
    mutate(count = n()) %>%
    dplyr::select(-SampleName) %>%
    dplyr::distinct(Group, count)
  
  missval <- overlaps %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::arrange(Group) %>%
    dplyr::select(-SampleName) %>%
    group_by(Protein, Group) %>% 
    summarise(Intensity = sum(Intensity)) %>% 
    left_join(annot_group, by="Group") %>% 
    mutate(status = 1-(Intensity/count)) %>%
    dplyr::select(Protein,status) %>%
    mutate(check = ifelse(status <= batch_missingness, 0, 1)) %>%
    group_by(Protein) %>% 
    summarise(check = sum(check))
  
  dropped = length(which(missval$check!=0))
  
  rownames(missval) <- missval$Protein
  frags_to_keep_c = row.names(missval)[which(missval$check==0)]
  fil_df = fil_df[(fil_df$Protein %in% frags_to_keep_c),]
  
  ##### LOGGING ######
  
  print("Drops from c)")
  print("Fragments dropped: ")
  print(dropped)
  print("Remaining Samples:")
  print(nrow(annotation))
  print("Remaining Fragments:")
  print(nrow(fil_df))
  print("\n")
  
  list_dfs <- list(fil_df, annotation)
  return(list_dfs)
}

#FullRunName in annotation file must be in the same order as input data
order_samples <- function(df, annotation){
  df = df[ , !(names(df) %in% c('ProteinName','PeptideSequence','FragmentIon','RT','Protein'))]
  sample_order = colnames(df)
  annotation <- annotation %>% arrange(factor(FullRunName, levels = sample_order))
  return(annotation)
}

# functions to plot normalized and unnorm box plot
calc_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(stats)
}

plot_boxnorm <- function(box_norm.m, color_list, exp_grp) {
  # necessary packages / libraries
  require(ggplot2)
  require(RColorBrewer)
  require(ggpubr)
  
  return(
    ggplot(data=box_norm.m, mapping=aes(x=FullRunName, y=Intensity, fill=get(exp_grp))) +
      #geom_boxplot(color="grey30") +
      stat_summary(fun.data = calc_stat, geom="boxplot") +
      theme_classic2()+
      scale_fill_manual(name=as.character(exp_grp), values = color_list[[exp_grp]]) +
      theme(text = element_text(size = 14,face="bold"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12,face="bold")) +
      theme(legend.position = "right") +
      labs(fill = exp_grp)
  )
}

# impute
nafunctions <- function(sample_matrix, method) {
  df <- df1 <- as.data.frame(sample_matrix)
  method <- tolower(method)
  if(method=="halfminimum"){
    for(i in 1:nrow(df)) {
      df[i,which(is.na(df[i,]))] <- min(df[i,], na.rm = T) / 2
    }
  }
  else if (method == "ranger") {
    df <- missRanger(df1, pmm.k = 3, splitrule = "extratrees", num.trees = 10, sample.fraction = 0.1, maxiter=10, seed=347)
  }
  else{
    stop("Unspported method")
  }
  df <- as.data.frame(df)
  return(df)
}

# function to create the df with annotations for the hca
hca_anno_df <- function(annotation_data, selected_annotations) {
  # generate df with annotations for heatmap
  df_with_annotations <- data.frame(matrix(ncol = 0, nrow = nrow(annotation_data)))
  
  for (i in selected_annotations) {
    print(i)
    annot_1 <- annotation_data[order(annotation_data[[i]]),][[i]]
    annot_1 <- list(i = annot_1)
    #print(annot_1)
    df_with_annotations[ , ncol(df_with_annotations) + 1] <- annot_1
    colnames(df_with_annotations)[ncol(df_with_annotations)] <- paste0(as.character(i))
  }
  return(df_with_annotations)
}

# heatmap function
plotly_heatmap <- function(for_heatmap, df_with_annotations, color_list) {
  # necessary libraries / packages
  require(plotly)
  require(iheatmapr)
  require(RColorBrewer)
  require(scales)
  require(stringr)
  
  set.seed(10) # consistent k-means
  
  p <- main_heatmap(
    for_heatmap,
    name = "Scaled<br>Intensity",
    colors = "RdYlBu",
    tooltip = setup_tooltip_options(
      value = TRUE,
      prepend_row = "Protein: ",
      prepend_col = "Sample: "
    ),
    layout = list(font = list(size = 0.5), height = 500)
  ) %>% 
    add_col_clustering(
      name = "<i>Sample</i><br>Cluster Assignment",
      method = "hclust",
      show_colorbar = F
    ) %>%
    add_col_annotation(
      colors = color_list,
      side = "bottom",
      df_with_annotations
    )
  p <- p %>% 
    add_col_title("Sample clustering",
                  side = "top",
                  size = 0.01,
                  buffer = 0,
                  font = list(size = 14)
    ) %>%
    # convert IHeatmapR object to Plotly object, needed to embed in HTML
    to_plotly_list() %>% 
    plotly::as_widget() %>% 
    config(
      displayModeBar = T,
      displaylogo = F
    ) %>%
    layout(height = 500, width = 800)
  return(p)
}

plotly_corr <- function(for_heatmap, anno_for_corr, col_element, replicate_filenames, color_list) {
  # necessary libraries / packages
  require(plotly)
  require(iheatmapr)
  require(RColorBrewer)
  require(scales)
  require(stringr)
  
  # sample order, separate for rows and columns
  sample_order = as.vector(anno_for_corr$FullRunName)
  sample_order2 = rev(sample_order)
  
  # annotation order, separate for rows and columns
  anno_order = as.vector(anno_for_corr[[col_element]])
  anno_df <- data.frame(matrix(ncol = 0, nrow = nrow(anno_for_corr)))
  anno_df[ , ncol(anno_df) + 1] <- anno_order
  colnames(anno_df)[ncol(anno_df)] <- paste0(as.character(col_element))
  
  anno_order2 = rev(anno_order)
  anno_df2 <- data.frame(matrix(ncol = 0, nrow = nrow(anno_for_corr)))
  anno_df2[ , ncol(anno_df2) + 1] <- anno_order2
  colnames(anno_df2)[ncol(anno_df2)] <- paste0(as.character(col_element))
  
  # make correlation matrix
  sample_matrix.corr <- for_heatmap[, which((names(for_heatmap) %in% replicate_filenames)==TRUE)]
  
  sample_matrix.corr <- cor(sample_matrix.corr, method = "spearman")
  sample_matrix.corr[is.na(sample_matrix.corr)] <- 0
  sample_matrix.corr[!is.finite(sample_matrix.corr)] <- 0
  
  # order correlation matrix by sample order for rows and columns
  sample_matrix.corr <- sample_matrix.corr[, sample_order2]
  
  sample_matrix.corr <- as.data.frame(sample_matrix.corr) %>%
    rownames_to_column(var="sample") %>%
    slice(match(sample_order, sample)) %>%
    column_to_rownames(var="sample")
  
  sample_matrix.corr <- as.matrix(sample_matrix.corr)
  
  # get corresponding color_list for correlation annotation
  # because sometimes we will not want to print all batches for correlation
  corr_anno_color_list <- c()
  for (i in seq(1,length(color_list[[col_element]]))){
    color <- color_list[[col_element]][i]
    #print(paste(names(color_list[[col_element]][i]), "goes", color))
    #print(color)
    if (names(color_list[[col_element]][i]) %in% anno_for_corr[[col_element]]) {
      corr_anno_color_list <- append(corr_anno_color_list, color)
    }
  }
  #corr_anno_color_list <- as.vector(corr_anno_color_list)
  #print(corr_anno_color_list)
  corr_anno_2 <- list()
  corr_anno_2[[col_element]] <- corr_anno_color_list
  #print(corr_anno_2)
  
  p <- main_heatmap(
    sample_matrix.corr,
    name = "Scaled<br>Intensity",
    colors = "RdYlBu",
    tooltip = setup_tooltip_options(
      value = FALSE,
      prepend_row = "FullRunName: "
    ),
    layout = list(font = list(size = 0.5), height = 500)
  ) %>% 
    add_col_annotation(
      colors = corr_anno_2,
      side = "top",
      anno_df
    ) %>%
    add_row_annotation(
      colors = corr_anno_2,
      side = "left",
      anno_df2
    )
  p <- p %>% 
    # convert IHeatmapR object to Plotly object, needed to embed in HTML
    to_plotly_list() %>% 
    plotly::as_widget() %>% 
    config(
      displayModeBar = T,
      displaylogo = F
    ) %>%
    layout(height = 500, width = 500)
  return(p)
}

################## function to add multiple linebreaks 
linebreaks <- function(n){HTML(strrep(br(), n))}

################### list of impute methods
methods <- c("halfminimum", "ranger")

################### ui
ui <- tagList(
  fluidPage(theme = shinytheme("sandstone"),
    shinyjs::useShinyjs(),
    tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
      }
    ")),
    tags$script(
      'function checkifrunning() {
        var is_running = $("html").attr("class").includes("shiny-busy");
        if (is_running){
          $("#loading").show()
        } else {
          $("#loading").hide()
        }
      }; 
      setInterval(checkifrunning, 50)'
    ), 
    tags$style(
      " body { text-align:left; }
      #loading {
        display: inline-block;
        border: 3px solid #f3f3f3; 
        border-top: 3px solid #3498db; 
        border-radius: 50%;
        width: 50px;
        height: 50px;
        animation: spin 1s ease-in-out infinite;
      }

      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }"
    ),
    tags$head(tags$script(type="text/javascript", src = "logo.js")),
    navbarPage(title="",
      id="maintab",
      tabPanel(
        "Home",
        uiOutput("homeui"),
      ), 
      tabPanel(
        "Settings",
        sidebarPanel(
          h4(
            "Upload Data + Select parameters",
            tags$span(
              id = 'span1',
              `data-toggle` = "tooltip",
              title = 'In this part, users can upload their own data following the instructions in the main panel. The example data can be automatically uploaded when users click "Load example data" below. The other parameters also get chosen with defaults using this example data.',
              tags$span(class = "glyphicon glyphicon-question-sign")
            )
          ),
          radioButtons(
            "data_type",
            label = "",
            choices = list("Load experimental data" = 1,"Load example data" = 2),
            selected = 1
          ),
          fileInput("unnorm_file", "1. Choose unnormalized file (tab delimited .txt) *", accept = ".txt"),
          span(uiOutput("unnorm_check"), style="color:red"),
          fileInput("norm_file", "2. Choose normalized file (tab delimited .txt)", accept = ".txt"),
          span(uiOutput("norm_check"), style="color:red"),
          selectInput("protein", "3. Choose the column containing the protein name *", choices=c(), multiple=FALSE, selected=NULL),
          fileInput("anno_file", "4. Choose annotation file (tab delimited .txt) *", accept = ".txt"),
          span(uiOutput("anno_check"), style="color:red"),
          selectInput("cols_of_int", "5. Choose the columns to correct for batch-effect *", choices=c(), multiple=TRUE, selected=NULL),
          selectInput("exp_grp", "6. Choose the column containing the experimental group *", choices=c(), multiple=FALSE, selected=NULL),
          span(uiOutput("exp_grp_check"), style="color:red"),
          selectInput("sample_names", "7. Choose the column containing the sample names (note that this column should match the headers in normalized data file) *", choices=c(), multiple=FALSE, selected=NULL),
          span(uiOutput("sample_names_check"), style="color:red"),
          actionButton("nextId", "Next", class = "btn-warning"),
        width = 3),
        mainPanel(
          h3("Input parameters"),
          p("On uploading your files, a preview will be available for review.", tags$b("The normalized and unnormalized file must contain a column with protein names and the remaining columns should be intensities corresponding to samples"), "as listed in the annotation file.", tags$b("The anotation file should contain columns to correct for batch-effect, a biological experimental group in which you want to retain variation, and a column with sample names that match the normalized/unnormalized file."), "Each of these columns can be specified in the menu on the left panel after uploading the data files."),
          linebreaks(1),
          
          h5(strong("Unnormalized file preview:")),
          dataTableOutput("unnorm_table"),
          linebreaks(2),
          hr(style = "border-top: 1px solid #000000;"),
          
          h5(strong("Normalized file preview:")),
          dataTableOutput("norm_table"),
          linebreaks(2),
          hr(style = "border-top: 1px solid #000000;"),
          
          h5(strong("Annotation file preview:")),
          dataTableOutput("anno_table"),
        ),
      ),
      tabPanel(
        "Initial analysis",
        navlistPanel(
          tabPanel("Sample distribution", 
                    h5(strong("Plot showing distribution of samples across batches")),
                    plotOutput("initial_plot1", width = "100%")  %>% withSpinner(color="#0dc5c1")), 
          
          tabPanel("Sample matrix",
                   h5(strong("Balloon plots showing how samples in experimental groups are distributed among the different batches")),
                   uiOutput("balloon_plots")  %>% withSpinner(color="#0dc5c1")),
          
          tabPanel("Missingness", 
                    h5(strong("Plot showing missingness distribution in overall data")),
                    plotOutput("pareto_plot")  %>% withSpinner(color="#0dc5c1"),
                    hr(style = "border-top: 1px solid #000000;"),
                   
                    h5(strong("Plots showing missingness distribution across batches chosen to correct on")),
                    uiOutput("missing_plots2")  %>% withSpinner(color="#0dc5c1"),
                    hr(style = "border-top: 1px solid #000000;"),
                   
                    h5(strong("Plot showing missingness distribution across experimental group")),
                    plotOutput("missing_plot3")  %>% withSpinner(color="#0dc5c1")),
          
          tabPanel("Take aways",
                    h5(strong("Stats based on sample distribution and missingness:")),
                    tableOutput("init_table")  %>% withSpinner(color="#0dc5c1"),
                    linebreaks(1),
                   
                    h5(strong("Importance of distribution of samples across batches:")),
                    htmlOutput("distri_text"),
                    linebreaks(1),
                   
                    h5(strong("Importance of missing data across samples and batches:")),
                    htmlOutput("missing_text"),
                    linebreaks(1),
                   
                    textOutput("missing_note1"),
                    textOutput("missing_note2"),
                    hr(style = "border-top: 1px solid #000000;"),
                   
                    h5(strong("Stats based on sample distribution per experimental group:")),
                    tableOutput("balloon_table")  %>% withSpinner(color="#0dc5c1"),
                    linebreaks(1),
                   
                    h5(strong("Importance of sample distribution within experimental group:")),
                    htmlOutput("balloon_text"),
                    hr(style = "border-top: 1px solid #000000;"),
                   
                    h5(strong("What's next:")),
                    textOutput("take_away_text"), 
                    linebreaks(1)),
        ),
      ),
      tabPanel(
        "Diagnosis+Filtering",
        sidebarPanel(
          radioButtons("impute_method", "Choose the method to be used for imputation*", methods),
          textInput("output_prefix", label = "Pick a prefix for output files*", value = "test"),
          textInput("iRT_prot", label = "Enter the name of your iRT protein as listed in your data file"),
          selectInput("samps_for_corr", "Choose the \"experimental group:batch\" combo for correlation. Pro tip: Use replicates here.", choices=c(), multiple=TRUE, selected=NULL),
          selectInput("var_to_correct_on", "Choose the variable you would like to correct on*:", choices=c()),
          sliderInput(inputId = "sample_threshold", 
                      label = "Sample Threshold*", 
                      value = 0.7, min = 0, max = 1),
          sliderInput(inputId = "expgrp_threshold", 
                      label = "Experimental group Threshold*", 
                      value = 0.5, min = 0, max = 1),
          sliderInput(inputId = "batch_threshold", 
                      label = "Batch Threshold*", 
                      value = 0.7, min = 0, max = 1),
          tags$div(style="display:inline-block",title="CAUTION: Once the submit button is hit the imputation method cannot be changed since re-imputing will take a lot of time. Be sure to select the appropriate option before hitting Submit.",actionButton("submitId2", "Submit"), class = "btn-warning"),
          #actionButton("submitId2", "Submit"),
          tags$div(style="display:inline-block",title="CAUTION: Once the continue button is hit the input parameters cannot be changed since filtering and batch-correction will take time. Be sure to provide the appropriate options before hitting Continue.",actionButton("continue", "Continue"), class = "btn-warning"),
          #actionButton("continue", "Continue"),
          actionButton("results", "Go to Results", class = "btn-warning"),
          width = 3),
        mainPanel(
          h3("Diagnosis for batch to correct on"),
          textOutput("results_msg"),
          h5(strong("PVCA before correction and filtering")),
          plotOutput("pvca_before_bc_init") %>% withSpinner(color="#0dc5c1"),
          linebreaks(2),
          textOutput("pvca_table"),
          hidden(hr(id = "line_pvca", style = "border-top: 1px solid #000000;")),
          
          hidden(h3(id = "heading_filt", "Diagnosis for filtering")),
          hidden(h5(id = "heading_perc_missing_plot", strong("Plot showing percent of missingness within each sample"))),
          plotOutput("filt_plot1")  %>% withSpinner(color="#0dc5c1"),
          hidden(hr(id = "line_filt_plot1", style = "border-top: 1px solid #000000;")),
          
          hidden(h5(id = "heading_missing_expgrp", strong("Missingness percent by experimental group"))),
          plotOutput("filt_plot2")  %>% withSpinner(color="#0dc5c1"),
          hidden(hr(id = "line_filt_plot2", style = "border-top: 1px solid #000000;")),
          
          hidden(h5(id = "heading_missing_batch", strong("Missingness percent by batch to correct on"))),
          plotOutput("filt_plot3")  %>% withSpinner(color="#0dc5c1"),
        ),
      ),
      tabPanel(
        "Results",
        navlistPanel(id="ResultsNav",
          tabPanel("Filtered results",
                   tabsetPanel(
                      tabPanel("Filtered stats",
                        h5(strong("Number of samples:")),
                        textOutput("init_samps"),
                        textOutput("filt_samps"),
                        hr(style = "border-top: 1px solid #000000;"),
                   
                        h5(strong("Number of fragments:")),
                        textOutput("init_frags"),
                        textOutput("filt_frags"),
                        hr(style = "border-top: 1px solid #000000;"),
                   
                        h5(strong("Missingness before/after filtering:")),
                        textOutput("filt_stats_note1"),
                        textOutput("filt_stats_note2"),
                        hr(style = "border-top: 1px solid #000000;")),
          
                      tabPanel("Sample distribution", 
                        h5(strong("Plot showing distribution of samples across batches (after filtering)")),
                        plotOutput("filtered_plot1", width = "100%")  %>% withSpinner(color="#0dc5c1")), 
          
                      tabPanel("Sample matrix",
                        h5(strong("Balloon plots showing how samples in experimental groups are distributed among the different batches (after filtering)")),
                        uiOutput("filtered_balloon_plots")  %>% withSpinner(color="#0dc5c1")),
                      
                      tabPanel("Missingness", 
                        h5(strong("Plots showing missingness distribution in overall data (after filtering)")),
                        plotOutput("filtered_plot2") %>% withSpinner(color="#0dc5c1"),
                        hr(style = "border-top: 1px solid #000000;"),
                   
                        h5(strong("Plots showing missingness distribution across batches chosen to correct on (after filtering)")),
                        uiOutput("filtered_plots3") %>% withSpinner(color="#0dc5c1"),
                        hr(style = "border-top: 1px solid #000000;"),
                   
                        h5(strong("Plot showing missingness distribution across experimental group (after filtering)")),
                        plotOutput("filtered_plot4") %>% withSpinner(color="#0dc5c1")),
                   ),
                  ),
          tabPanel("Normalization",
                   tabsetPanel(
                      tabPanel("Unnormalized", plotOutput("unnorm_box") %>% withSpinner(color="#0dc5c1")),
                      tabPanel("Normalized", plotOutput("norm_box") %>% withSpinner(color="#0dc5c1")),
                   ),
                  ),
          tabPanel("PVCA",
                   tabsetPanel(
                     tabPanel("Before correction", textOutput("pvca_before_bc_note"), plotOutput("pvca_before_bc") %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat only correction", plotOutput("pvca_after_across_bc") %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat and loess correction", plotOutput("pvca_after_bc") %>% withSpinner(color="#0dc5c1"))
                   )
                  ), 
          tabPanel("PCA",
                   tabsetPanel(
                      tabPanel("Before correction", plotOutput("pca_before_bc") %>% withSpinner(color="#0dc5c1")), 
                      tabPanel("After combat only correction", plotOutput("pca_after_across_bc") %>% withSpinner(color="#0dc5c1")), 
                      tabPanel("After combat and loess correction", plotOutput("pca_after_bc") %>% withSpinner(color="#0dc5c1"))
                   )
                  ),
          #tabPanel("HCA",
          #         tabsetPanel(
          #           tabPanel("Before correction", plotlyOutput("hca_before_bc") %>% withSpinner(color="#0dc5c1")), 
          #           tabPanel("After combat only correction", plotlyOutput("hca_after_across_bc") %>% withSpinner(color="#0dc5c1")), 
          #           tabPanel("After combat and loess correction", plotlyOutput("hca_after_bc") %>% withSpinner(color="#0dc5c1"))
          #         )
          #        ),
          tabPanel("iRT mapping",
                   tabsetPanel(id="iRT mapping",
                     tabPanel("Before correction", plotOutput("irt_before_bc") %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat only correction", plotOutput("irt_after_across_bc") %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat and loess correction", plotOutput("irt_after_bc") %>% withSpinner(color="#0dc5c1"))
                   )
                  ),
          tabPanel("Correlation",
                   tabsetPanel(id="Correlation",
                     tabPanel("Before correction", plotlyOutput("corr_before_bc", height = 500, width = 700) %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat only correction", plotlyOutput("corr_after_across_bc", height = 500, width = 700) %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat and loess correction", plotlyOutput("corr_after_bc", height = 500, width = 700) %>% withSpinner(color="#0dc5c1"))
                   )
                  ),
          tabPanel("Downloads",
                   h5(strong("Report:")),
                   downloadButton("report", "Download report") %>% withSpinner(color="#0dc5c1"),
                   linebreaks(2),
                   hr(style = "border-top: 1px solid #000000;"),
                   
                   h5(strong("Primary files for down-stream analysis:")),
                   textOutput("bc_type_selection_note"),
                   linebreaks(2),
                   downloadButton("after_preimpute", "After complete batch-correction, pre-imputed file"),
                   linebreaks(2),
                   downloadButton("after_postimpute", "After complete batch-correction, post-imputation file"), 
                   linebreaks(2),
                   downloadButton("after_combat_preimpute", "After combat batch-correction, pre-imputed file"),
                   linebreaks(2),
                   downloadButton("after_combat_postimpute", "After combat batch-correction, post-imputation file"),
                   linebreaks(2),
                   hr(style = "border-top: 1px solid #000000;"),
                   
                   
                   h5(strong("Secondary files:")),
                   downloadButton("filtnorm", "Filtered normalized file"),
                   linebreaks(2),
                   downloadButton("filtunnorm", "Filtered unnormalized file"),
                   linebreaks(2),
                   downloadButton("filtanno", "Filtered annotation file"),
                   linebreaks(2),
                   downloadButton("before_preimpute", "Before batch-correction, pre-imputed file"),
                   linebreaks(2),
                   downloadButton("before_postimpute", "Before batch-correction, post-imputation file"))
        )
      )
    )
  )
)

################### server
server <- function(input, output, session) {
  # increase file upload size limit to 60 MB
  options(shiny.maxRequestSize=90*1024^2)
  
  # To validate input
  iv <- InputValidator$new()
  iv$add_rule("irt_prot", sv_optional())
  iv$add_rule("samps_for_corr", sv_optional())
  iv$add_rule(
    "output_prefix",
    sv_regex("^[a-zA-Z0-9]*$", "Only alphanumeric characters allowed")
  )
  
  # `enable()` the validation rules and 
  iv$enable()
  
  # disable tabs
  shinyjs::disable(selector = '.navbar-nav a[data-value="Initial analysis"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Diagnosis+Filtering"')
  shinyjs::disable(selector = '.navbar-nav a[data-value="Results"')
  
  # disable file upload buttons if example data is chosen 
  observe({
    if (input$data_type == 2) {
      shinyjs::disable("norm_file")
      shinyjs::disable("unnorm_file")
      shinyjs::disable("anno_file")
      shinyjs::disable("protein")
      shinyjs::disable("anno_file")
      shinyjs::disable("cols_of_int")
      shinyjs::disable("exp_grp")
      shinyjs::disable("sample_names")
    }
    else {
      shinyjs::enable("norm_file")
      shinyjs::enable("unnorm_file")
      shinyjs::enable("anno_file")
      shinyjs::enable("protein")
      shinyjs::enable("anno_file")
      shinyjs::enable("cols_of_int")
      shinyjs::enable("exp_grp")
      shinyjs::enable("sample_names")
    }
  })
  
  ########### home 
  output$homeui <- renderUI({
    fluidRow(
      div(
        id="mainbody",
        column(2),
        column(
          8,
          div(style="text-align:center;margin-top:20px;font-size:140%;color:darkorange",
              HTML("<em>Welcome to BIRCH</em>")),
          div(style="width:fit-content;width:-webkit-fit-content;width:-moz-fit-content;font-size:100%;margin-top:10px",
              HTML("<b>BIRCH</b> (Batch-effect Identification, Rectification, and Conclusive-anlaysis on Heterogenous data) is a web-app that can be used for reducing batch-effect in proteomics data. It generally becomes necessary to correct for batch-effect in large datasets since processing steps such as sample preparation and data acquisition tends to add noise to the data, that in-turn effects biological conclusions. This tool aims at keeping meaningful biological variation, while simultaneously reducing batch-effect due to other external factors. The steps involved in individual tabs are as follows:")),
          div(style="text-align:center;margin-top: 20px",
              a(href='#',
                img(src = "workflow.PNG", height = 300, width = 600))),
          div(style="text-align:center;margin-top: 20px",
              a(href="https://github.com/csmc-vaneykjlab/BatchCorrectionTool", class="btn btn-default",
                "GitHub"),
              a(href="https://github.com/csmc-vaneykjlab/BatchCorrectionTool", class="btn btn-default",
                "Cite us"),
              a(href="mailto:ArchanaSubrama.Bhat@cshs.org", class="btn btn-default",
                "Email")),
          div(style="text-align:center;margin-top: 20px",
              HTML("")),
          ),
        column(2)
        )
      )
  })
  
  
  ########### Check unnorm file and display it
  output$unnorm_table <- renderDataTable({
    if (input$data_type == 2) {
      # Read example data and display it
      data_unnorm <- read.delim("example_unnorm.txt", header = TRUE)
      return(datatable(data_unnorm, options = list(pageLength = 10)))
    }
    else {
      # validate issues in file
      validate(
        need(!is.null(input$unnorm_file),
             "Please upload unnormalized data file")
      )
      
      # read file once provided
      file <- input$unnorm_file
      data_unnorm <- read.delim(file$datapath, header = TRUE)
      check <- select_if(data_unnorm, is.numeric)
      ncol_unnorm <- ncol(data_unnorm)
      ncol_check <- ncol(check)
      
      # validate issues in file
      # make sure only one text column exists
      validate(
        need(ncol_check == (ncol_unnorm-1),
             "Please upload unnormalized data file with only one protein column (containing text), and the rest as intensities (numeric columns) for each sample!")
      )
      
      # If everything passes, print file
      return(datatable(data_unnorm, options = list(pageLength = 10))) 
    }
  })
  
  output$unnorm_check <- renderUI({
    req(input$unnorm_file)
    
    # read file once provided
    file <- input$unnorm_file
    data_unnorm <- read.delim(file$datapath, header = TRUE)
    check <- select_if(data_unnorm, is.numeric)
    ncol_unnorm <- ncol(data_unnorm)
    ncol_check <- ncol(check)
    
    if(ncol_check != (ncol_unnorm-1)) {
      return(p("Please upload unnormalized data file with only one protein column (containing text), and the rest as intensities (numeric columns) for each sample!"))
    } else {
      return("")
    }
  })
  
  unnorm_file_name <- reactive({ 
    file <- input$unnorm_file
    return(file$datapath)
  })
  
  ############# Process unnorm file
  sample_matrix_unnorm_init <- reactive({
    if (input$data_type == 2) {
      # Read example data and display it
      unnorm_dataframe <- read.delim("example_unnorm.txt", header = TRUE)
      return(unnorm_dataframe)
    }
    else {
      req(input$unnorm_file)
      unnorm_dataframe <- input$unnorm_file
      unnorm_dataframe <- read.delim(unnorm_dataframe$datapath, header = TRUE)
      return(unnorm_dataframe) 
    }
  })
  
  ########### Check norm file and display it
  output$norm_table <- renderDataTable({
    if (input$data_type == 2) {
      # Read example data and display it
      data_norm <- read.delim("example_norm.txt", header = TRUE)
      return(datatable(data_norm, options = list(pageLength = 10)))
    }
    else {
      # validate issues in file
      validate(
        need(!is.null(input$norm_file),
             "Please upload normalized data file. This is an optional argument. If you do not have a normalized file, quantile normalization will be performed and the normalized data will be used for further analysis. You can download the normalized file from the Results section.")
      )
      
      # using norm file again
      file_norm <- input$norm_file
      data_norm <- read.delim(file_norm$datapath, header = TRUE)
      check_norm <- select_if(data_norm, is.numeric)
      ncol_norm <- ncol(data_norm)
      ncol_check_norm <- ncol(check_norm)
      nrow_norm <- nrow(data_norm)
      
      # read norm file again to make sure ncol, nrow in norm and unnorm are same
      data_unnorm <- sample_matrix_unnorm_init()
      check_unnorm <- select_if(data_unnorm, is.numeric)
      ncol_unnorm <- ncol(data_unnorm)
      ncol_check_unnorm <- ncol(check_unnorm)
      nrow_unnorm <- nrow(data_unnorm)
      
      # validate issues in file
      validate(
        need(ncol_check_norm == (ncol_norm-1),
             "Please upload normalized data file with only one protein column (containing text), and the rest as intensities (numeric columns) for each sample!") %then%
          need(ncol_norm == ncol_unnorm,
               "Make sure number of columns/samples are same in unnormalized and normalized files, kindly re-load the data accordingly.") %then%
          need(nrow_norm == nrow_unnorm,
               "Make sure number of rows/proteins are same in unnormalized and normalized files, kindly re-load the data accordingly.")
      )
      
      # If everything passes, print file
      return(datatable(data_norm, options = list(pageLength = 10))) 
    }
  })
  
  output$norm_check <- renderUI({
    req(input$norm_file)
    
    # using norm file again
    file_norm <- input$norm_file
    data_norm <- read.delim(file_norm$datapath, header = TRUE)
    check_norm <- select_if(data_norm, is.numeric)
    ncol_norm <- ncol(data_norm)
    ncol_check_norm <- ncol(check_norm)
    nrow_norm <- nrow(data_norm)
    
    # read norm file again to make sure ncol, nrow in norm and unnorm are same
    data_unnorm <- sample_matrix_unnorm_init()
    check_unnorm <- select_if(data_unnorm, is.numeric)
    ncol_unnorm <- ncol(data_unnorm)
    ncol_check_unnorm <- ncol(check_unnorm)
    nrow_unnorm <- nrow(data_unnorm)
    
    if(ncol_check_norm != (ncol_norm-1)) {
      return(p("Please upload normalized data file with only one protein column (containing text), and the rest as intensities (numeric columns) for each sample!"))
    } else if (ncol_norm != ncol_unnorm) {
      return(p("Make sure number of columns and sample names are same in unnormalized and normalized files, kindly re-load the data accordingly."))
    } else if (nrow_norm != nrow_unnorm) {
      return(p("Make sure number of rows and protein names are same in unnormalized and normalized files, kindly re-load the data accordingly."))
    }else {
      return("")
    }
  })
  
  norm_file_name <- reactive({ 
    file <- input$norm_file
    return(file$datapath)
  })
  
  ############### Process norm file
  sample_matrix_init <- reactive({
    if (input$data_type == 2) {
      # Read example data and display it
      norm_dataframe <- read.delim("example_norm.txt", header = TRUE)
      return(norm_dataframe)
    }
    else {
      req(input$norm_file)
      norm_dataframe <- input$norm_file
      norm_dataframe <- read.delim(norm_dataframe$datapath, header = TRUE)
      return(norm_dataframe) 
    }
  })
  
  ########### Display anno file
  output$anno_table <- renderDataTable({
    if (input$data_type == 2) {
      # Read example data and display it
      data_anno <- read.delim("example_anno.txt", header = TRUE)
      return(datatable(data_anno, options = list(pageLength = 10)))
    }
    else {
      validate(
        need(!is.null(input$anno_file),
             "Please upload annotation file")
      )
      
      # read norm file again to make sure ncol, nrow in norm and unnorm are same
      data_unnorm <- sample_matrix_unnorm_init()
      ncol_unnorm <- ncol(data_unnorm)
      
      # read anno file
      file_anno <- input$anno_file
      data_anno <- read.delim(file_anno$datapath, header = TRUE)
      nrow_anno <- nrow(data_anno)
      
      # validate issues in file
      validate(
        need((ncol_unnorm-1) == nrow_anno,
             "Please ensure that the number of samples are same in annotation (rows) and normalized/unnormalized protein intensities file (columns)! Re-load your files as necessary.")
      )
      
      # If everything passes, print file
      return(datatable(data_anno, options = list(pageLength = 10))) 
    }
  })
  
  output$anno_check <- renderUI({
    req(input$anno_file)
    
    # read norm file again to make sure ncol, nrow in norm and unnorm are same
    data_unnorm <- sample_matrix_unnorm_init()
    ncol_unnorm <- ncol(data_unnorm)
    
    # read anno file
    file_anno <- input$anno_file
    data_anno <- read.delim(file_anno$datapath, header = TRUE)
    nrow_anno <- nrow(data_anno)
    
    if((ncol_unnorm-1) != nrow_anno) {
      return(p("Please ensure that the number of samples are same in annotation (rows) and normalized/unnormalized protein intensities file (columns)! Re-load your files as necessary."))
    } else {
      return("")
    }
  })
  
  anno_file_name <- reactive({ 
    file <- input$anno_file
    return(file$datapath)
  })
  
  ########### Process anno file
  annotation_data_init <- reactive({
    if (input$data_type == 2) {
      # Read example data and display it
      anno_dataframe <- read.delim("example_anno.txt", header = TRUE)
      return(anno_dataframe)
    }
    else {
      req(input$anno_file)
      anno_dataframe <- input$anno_file
      anno_dataframe <- read.delim(anno_dataframe$datapath, header = TRUE)
      return(anno_dataframe) 
    }
  })
  
  ########### Process cols_of_int
  observe({
    if (input$data_type == 2) {
      updateSelectInput(session, "cols_of_int",
                        choices = colnames(annotation_data_init()),
                        selected = "Digestion_batch") 
    }
    else {
      updateSelectInput(session, "cols_of_int",
                        choices = colnames(annotation_data_init()),
                        selected = colnames(annotation_data_init()[0])) 
    }
  })
  
  rv <- reactiveVal(NULL)
  observeEvent(input$cols_of_int, {
    rv(input$cols_of_int)
  })
  
  ########## Process exp_grp
  observe({
    if (input$data_type == 2) {
      updateSelectInput(session, "exp_grp",
                        choices = colnames(annotation_data_init()),
                        selected = "ExperimentalGroup")
    }
    else {
      updateSelectInput(session, "exp_grp",
                        choices = colnames(annotation_data_init()),
                        selected = "")
      #selected = colnames(annotation_data_init()[0])) 
    }
  })
  
  exp_grp_val <- reactiveVal(NULL)
  observeEvent(input$exp_grp, {
    exp_grp_val(input$exp_grp)
  })
  
  # Print an error in side panel if exp grp is in cols of int
  output$exp_grp_check <- renderUI({
    req(input$anno_file, input$exp_grp, input$cols_of_int)
    
    exp_grp <- exp_grp_val()
    col_of_int_vec <- as.vector(rv())
    
    if(exp_grp %in% col_of_int_vec) {
      return(p("Please ensure that the experimental group is not one of the columns to correct on. Either pick another experimental group, or change your columns of interest."))
    } else {
      return("")
    }
  })
  
  ########## Process sample_names
  observe({
    if (input$data_type == 2) {
      updateSelectInput(session, "sample_names",
                        choices = colnames(annotation_data_init()),
                        selected = "SampleName")
    }
    else {
      updateSelectInput(session, "sample_names",
                        choices = colnames(annotation_data_init()),
                        selected = "")
      #selected = colnames(annotation_data_init()[0])) 
    }
  })
  
  sample_names_val <- reactiveVal(NULL)
  observeEvent(input$sample_names, {
    sample_names_val(input$sample_names)
  })
  
  ########## Process protein
  observe({
    if (input$data_type == 2) {
      updateSelectInput(session, "protein",
                        choices = colnames(sample_matrix_unnorm_init()),
                        selected = "ProteinName")
    }
    else {
      updateSelectInput(session, "protein",
                        choices = colnames(sample_matrix_unnorm_init()),
                        selected = "")
      #selected = colnames(sample_matrix_init()[0])) 
    }
  })
  
  protein_val <- reactiveVal(NULL)
  observeEvent(input$protein, {
    protein_val(input$protein)
  })
  
  output$sample_names_check <- renderUI({
    req(input$anno_file, input$unnorm_file, input$protein, input$sample_names)
    
    protein <- protein_val()
    unnorm_col_names <- as.vector(colnames(sample_matrix_unnorm_init()))
    unnorm_col_names <- unnorm_col_names[! unnorm_col_names %in% c(protein)]
    
    anno_data <- annotation_data_init()
    sample_names <- sample_names_val()
    
    sample_names_vec <- as.vector(anno_data[[sample_names]])
    sample_names_fixed <- c()
    for (each in sample_names_vec) {
      first_char <- substr(each, 1, 1)
      # if first character cannot be converted to integer then is.na will return TRUE and the first character is non-numeric in nature
      first_char_check <- is.na(as.integer(first_char))
      if (first_char_check==TRUE) {
        new_name <- each
      } else {
        new_name <- paste0('X', each)
      }
      sample_names_fixed <- append(sample_names_fixed, new_name)
    }
    
    if (all(unnorm_col_names==sample_names_fixed)) {
      return("")
    } else {
      return(p("Please ensure that the sample names in unnorm/norm file are same as the sample names in annotation file."))
    }
  })
  
  # activate the next button once all required input is provided
  observe({
    if (input$data_type == 2) {
      if (!isTruthy(input$protein) || !isTruthy(input$exp_grp) || is.null(input$cols_of_int) || !isTruthy(input$sample_names))  {
        shinyjs::disable("nextId")
      } else {
        shinyjs::enable("nextId")
      } 
    }
    else {
      if (is.null(input$unnorm_file) || is.null(input$anno_file) || !isTruthy(input$protein) || !isTruthy(input$exp_grp) || is.null(input$cols_of_int) || !isTruthy(input$sample_names))  {
        shinyjs::disable("nextId")
      } else {
        shinyjs::enable("nextId")
      } 
    }
  })
  
  ########### Process files - part 2 - calling functions to clean and check the data
  annotation_data <- reactive({
    req(input$sample_names)
    #req(input$anno_file, input$sample_names)
    anno_dataframe <- annotation_data_init()
    sample_names <- sample_names_val()
    
    # make sure the sample names starting with digits are prefixed with "X"
    sample_names_vec <- as.vector(anno_dataframe[[sample_names]])
    sample_names_fixed <- c()
    for (each in sample_names_vec) {
      first_char <- substr(each, 1, 1)
      # if first character cannot be converted to integer then is.na will return TRUE and the first character is non-numeric in nature
      first_char_check <- is.na(as.integer(first_char))
      if (first_char_check==TRUE) {
        new_name <- each
      } else {
        new_name <- paste0('X', each)
      }
      sample_names_fixed <- append(sample_names_fixed, new_name)
    }
    
    # replace the sample names with the fixed and prefixed sample names
    anno_dataframe[[sample_names]] <- sample_names_fixed
    
    # call methods to clean the anno file
    anno_dataframe <- colClean2(anno_dataframe, sample_names)
    anno_dataframe <- rename_level3(anno_dataframe, sample_names)
    
    # create a column called order if it doesn't exist
    col_names <- names(anno_dataframe)
    if ("order" %in% col_names) {
      return(anno_dataframe) 
    } else {
      anno_dataframe$order <- 1:nrow(anno_dataframe)
      return(anno_dataframe) 
    }
  })
  
  sample_matrix <- reactive({
    #req(input$unnorm_file)
    if (input$data_type == 2) {
      norm_dataframe <- sample_matrix_init()
      norm_dataframe <- colClean(norm_dataframe)
      return(norm_dataframe) 
    }
    else {
      if (is.null(input$norm_file)) {
        norm_dataframe <- sample_matrix_unnorm_init()
        norm_dataframe <- colClean(norm_dataframe)
        protein <- protein_val()
        rownames(norm_dataframe) <- norm_dataframe[[protein]]
        norm_dataframe <- norm_dataframe[,-which(colnames(norm_dataframe) %in% c(protein))] 
        # quantile normalization
        quantile_normalized_matrix <- normalize_data_dm(as.matrix(norm_dataframe), normalize_func = "quantile")
        mode(quantile_normalized_matrix) = "numeric"
        #quantile_normalized_matrix <- cbind(rownames(quantile_normalized_matrix), quantile_normalized_matrix) # unset rownames as protein names
        #colnames(quantile_normalized_matrix)[1] <- protein
        #rownames(quantile_normalized_matrix)<-NULL
        quantile_normalized_matrix <- as.data.frame(quantile_normalized_matrix)
        quantile_normalized_matrix <- tibble::rownames_to_column(quantile_normalized_matrix, protein)
        return(quantile_normalized_matrix)  
      } 
      else {
        norm_dataframe <- sample_matrix_init()
        norm_dataframe <- colClean(norm_dataframe)
        return(norm_dataframe) 
      } 
    }
  })

  sample_matrix_unnorm <- reactive({
    #req(input$unnorm_file)
    #unnorm_dataframe <- input$unnorm_file
    unnorm_dataframe <- sample_matrix_unnorm_init()
    unnorm_dataframe <- colClean(unnorm_dataframe)
    return(unnorm_dataframe) 
  })
  
  
  ###### Activate initial anal and filtering tab
  observeEvent(input$nextId, {
    updateNavbarPage(
      inputId = "maintab",
      selected = "Initial analysis"
    )
    shinyjs::enable(selector = '.navbar-nav a[data-value="Initial analysis"')
    shinyjs::enable(selector = '.navbar-nav a[data-value="Diagnosis+Filtering"')
  })
  
  ########### Generate initial analysis results
  input_plot1_reac <- reactive({
    req(input$cols_of_int, input$exp_grp, input$sample_names)
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
    technical_factors <- as.character(rv())
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors,technical_factors)
    print(selected_annotations)
    
    index = 0
    p_list <- list()
    for (i in selected_annotations) {
      index = index + 1
      print(i)
      print(index)
      #cat("### Samples by ", i, "\n\n", sep="")
      p_list[[index]] <- plot_bar(annotation_data(), i, color_list)
      #print(p)
    }
    
    return(grid.arrange(grobs=p_list, ncol=3))
  })
  
  output$initial_plot1 <- renderPlot({
    p <- input_plot1_reac()
    print(p)
  }, height = 400, width = 1000)

  
  # For missingness in samples
  # Make table with sample_matrix and annotation combined
  joined_df <- reactive({
    req(input$cols_of_int)
    #req(input$unnorm_file, input$cols_of_int)
    table2 <- sample_matrix() %>%
      select(everything()) %>%
      summarise_all(funs(sum(is.na(.)))) %>%
      rownames_to_column() %>%
      pivot_longer(cols=-rowname)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    
    joined_df <- merge(table2, annotation_data(), by.x = "name", by.y = "FullRunName")
    joined_df <- joined_df %>% 
      select(check_cols_of_interest, exp_grp_val(), "name", "value")
    
    return(joined_df)
  })
  
  # Balloon plot
  balloon_plots_reac <- reactive({
    balloon_plot_output_list <- lapply(1:length(as.character(rv())), function(i) {
      plotname <- paste("balloonplot", i, sep="")
      plotOutput(plotname, height=500)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items to display properly.
    do.call(tagList, balloon_plot_output_list)
  })
  
  num_plots <- 5
  for (i in 1:num_plots) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("balloonplot", my_i, sep="")
      output[[plotname]] <- renderPlot({
        req(input$cols_of_int)
        
        index <- as.list(rv()[my_i])
        index_char <- toString(index)
        print(index_char)
        
        exp_grp <- exp_grp_val()
        balloon_plot <- plot_balloon(annotation_data(), index_char, exp_grp)
        return(balloon_plot)
      })
    })
  }
  
  output$balloon_plots <- renderUI({
    p <- balloon_plots_reac()
    print(p)
  })
  
  balloon_plot_reac_report <- reactive({
    req(input$cols_of_int)

    balloon_plots <- c()
    
    technical_factors <- unlist(strsplit(as.character(rv()), " "))
    
    for (i in 1:length(technical_factors)) {
      
      balloon_plots[[i]] <- local({
        i <- i
        index_char <- technical_factors[[i]]
        print(index_char)
        
        exp_grp <- exp_grp_val()
        balloon_plot <- plot_balloon(annotation_data(), index_char, exp_grp)
        print(balloon_plot)
      })
    }
    return(balloon_plots)
  })
  
  # Pareto plot
  pareto_plot_reac <- reactive({
    # hard coding batch to first batch in cols of int 
    index <- as.list(rv()[1])
    index_char <- toString(index)
    protein <- protein_val()
    sample_matrix <- sample_matrix()
    overlaps <- call_missing_vals(sample_matrix, annotation_data(), index_char, protein)
    plot <- plot_pareto(overlaps, protein)
    return(plot)
  })
  
  output$pareto_plot <- renderPlot({
    p <- pareto_plot_reac()
    print(p)
  })
  
  # Heatmaps showing missingness
  # Insert the right number of plot output objects into the web page
  missing_plots2_reac <- reactive({
    plot_output_list2 <- lapply(1:length(as.character(rv())), function(i) {
      plotname <- paste("missplot", i, sep="")
      plotOutput(plotname)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items to display properly.
    do.call(tagList, plot_output_list2)
  })
  
  num_plots <- 5
  for (i in 1:num_plots) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("missplot", my_i, sep="")
      output[[plotname]] <- renderPlot({
        req(input$cols_of_int)
        
        index <- as.list(rv()[my_i])
        index_char <- toString(index)
        print(index_char)
        
        num_of_file <- length(rownames(sample_matrix()))
        if ((10/num_of_file)*40 > 10){
          FONTSIZE_sample_axis = 10
        }else{
          FONTSIZE_sample_axis = (10/num_of_file)*150
        }
        
        exp_grp <- exp_grp_val()
        sample_names <- sample_names_val()
        color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
        protein <- protein_val()
        overlaps <- call_missing_vals(sample_matrix(), annotation_data(), index_char, protein)
        miss_plot <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, index_char, protein)
        return(miss_plot)
      })
    })
  }
  
  output$missing_plots2 <- renderUI({
    p <- missing_plots2_reac()
    print(p)
  })
  
  missing_plots2_reac_report <- reactive({
    req(input$cols_of_int)
    
    missing_plots <- vector('list')
    
    technical_factors <- unlist(strsplit(as.character(rv()), " "))
    
    for (i in 1:length(technical_factors)) {
      
      missing_plots[[i]] <- local({
        i <- i
        index_char <- technical_factors[[i]]
        print(index_char)
        
        num_of_file <- length(rownames(sample_matrix()))
        if ((10/num_of_file)*40 > 10){
          FONTSIZE_sample_axis = 10
        }else{
          FONTSIZE_sample_axis = (10/num_of_file)*150
        }
        
        exp_grp <- exp_grp_val()
        sample_names <- sample_names_val()
        color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
        protein <- protein_val()
        overlaps <- call_missing_vals(sample_matrix(), annotation_data(), index_char, protein)
        miss_plot <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, index_char, protein)
      })
    }
    return(missing_plots)
  })
  
  missing_plot3_reac <- reactive({
    req(input$cols_of_int)
    
    num_of_file <- length(rownames(sample_matrix()))
    if ((10/num_of_file)*40 > 10){
      FONTSIZE_sample_axis = 10
    }else{
      FONTSIZE_sample_axis = (10/num_of_file)*150
    }
    
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix(), annotation_data(), exp_grp, protein)
    miss_plot2 <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, exp_grp, protein)
    
    return(miss_plot2)
  })
  
  output$missing_plot3 <- renderPlot({
    p <- missing_plot3_reac()
    print(p)
  })
  
  init_table_reac <- reactive({
    req(input$cols_of_int)
    #req(input$cols_of_int, input$anno_file)
    
    technical_factors <- unlist(strsplit(as.character(rv()), " "))
    init_table <- data.frame(matrix(ncol = 12, nrow=length(technical_factors)))
    x <- c("Batch_column", "Total_num_of_plates", "Max_num_of_samples", "Min_num_of_samples", "Variation_in_sample_distribution", "Plates_passing_distri_cutoff_init", "Plates_passing_distri_cutoff", "Max_missing_by_plate", "Min_missing_by_plate", "Variation_in_missingness", "Plates_passing_missing_cutoff_init", "Plates_passing_missing_cutoff")
    colnames(init_table) <- x
    init_table_index = 0
    
    protein <- protein_val()
    
    for (i in technical_factors) {
      init_table_index = init_table_index + 1  
      
      samps_table <- annotation_data() %>% 
        group_by("Group"=annotation_data()[[i]]) %>% 
        summarise(Sample_count = n()) %>%
        mutate(Sample_percentage=100*Sample_count/sum(Sample_count))
      
      overlaps_for_init_table <- pivot_longer(sample_matrix(), cols = c(everything(), -protein), names_to = "FullRunName", values_to = "Intensity") %>% 
        left_join(annotation_data(), by="FullRunName") %>% 
        dplyr::rename("SampleName" = "FullRunName" , "Group" := !!i) %>% 
        subset(select=c(protein,"Intensity","SampleName","Group")) %>%
        group_by(Group) %>%
        summarise(na_count = sum(is.na(Intensity)), not_na = sum(!is.na(Intensity))) %>%
        mutate(total = na_count + not_na) %>%
        mutate(perc_of_na = na_count/total*100) 
      
      overlaps_final <- merge(overlaps_for_init_table, samps_table)
      
      Batch_col <- toString(i)
      Max_missing_by_plate <- overlaps_final %>% summarise(max(perc_of_na))
      Min_missing_by_plate <- overlaps_final %>% summarise(min(perc_of_na))
      Variation_in_missingness <- (Max_missing_by_plate - Min_missing_by_plate)/Max_missing_by_plate
      Plates_passing_missing_cutoff <- overlaps_final %>% summarise(sum(perc_of_na < 50)) 
      Total_num_of_plates <- nrow(overlaps_final)
      Num_of_plates_passing_missing_cutoff <- paste0(Plates_passing_missing_cutoff, " out of ", Total_num_of_plates)
      
      Max_num_of_samples <- overlaps_final %>% summarise(max(Sample_count))
      Min_num_of_samples <- overlaps_final %>% summarise(min(Sample_count))
      Variation_in_sample_distribution <- (Max_num_of_samples - Min_num_of_samples)/Max_num_of_samples*100
      Plates_passing_distri_cutoff <- overlaps_final %>% summarise(sum(Sample_count > 25)) 
      Num_of_plates_passing_distri_cutoff <- paste0(Plates_passing_distri_cutoff, " out of ", Total_num_of_plates)
      
      row_to_add = c(Batch_col, Total_num_of_plates, Max_num_of_samples, Min_num_of_samples, Variation_in_sample_distribution, Plates_passing_distri_cutoff, Num_of_plates_passing_distri_cutoff, Max_missing_by_plate, Min_missing_by_plate, Variation_in_missingness, Plates_passing_missing_cutoff, Num_of_plates_passing_missing_cutoff)
      
      #init_table = rbind(init_table, row_to_add)
      init_table[init_table_index,] <- row_to_add
      #init_table[nrow(init_table)+1,] <- row_to_add
    }
    
    init_table <- init_table %>% 
      select(-Max_missing_by_plate, -Min_missing_by_plate, -Max_num_of_samples, -Min_num_of_samples, -Plates_passing_missing_cutoff_init, -Plates_passing_distri_cutoff_init, -Total_num_of_plates) %>%
      mutate(Variation_in_missingness = round(Variation_in_missingness, 4)) %>%
      mutate(Variation_in_sample_distribution = round(Variation_in_sample_distribution, 4))

    init_table$Variation_in_missingness <- paste0(init_table$Variation_in_missingness, " %")
    init_table$Variation_in_sample_distribution <- paste0(init_table$Variation_in_sample_distribution, " %")
    
    return(init_table)
  })
  
  output$init_table <- renderTable({
    init_table_reac()
  })
  
  balloon_table_reac <- reactive({
    req(input$exp_grp, input$cols_of_int)
    #req(input$exp_grp, input$cols_of_int, input$anno_file)
    exp_grp <- exp_grp_val()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    
    balloon_table <- data.frame(matrix(ncol = 3, nrow=length(check_cols_of_interest)))
    x <- c("Batch_column", "Plates_with_NO_samples", "Plates_with_ALL_samples")
    colnames(balloon_table) <- x
    balloon_table_index = 0
    
    for (batch_col in check_cols_of_interest) {
      balloon_table_index = balloon_table_index + 1
      anno_for_table <- annotation_data() %>% dplyr::count(across(exp_grp), across(batch_col), .drop = FALSE)
      
      anno_zeros <- anno_for_table %>%
        group_by(across(all_of(batch_col))) %>%
        summarize(count_zero = count(n==0))
      
      plates_with_no_samps_per_exp_grp <- sum(anno_zeros$count_zero > 0)
      
      anno_not_zeros <- anno_for_table %>%
        group_by(across(all_of(batch_col))) %>%
        summarize(count_not_zero = count(n!=0))
      
      plates_with_all_samps_per_exp_grp <- sum(anno_not_zeros$count_not_zero == 1)
      num_of_plates <- length(unique(annotation_data()[[batch_col]]))
      
      plates_with_no_samps_per_exp_grp_to_print <- paste0(plates_with_no_samps_per_exp_grp, " out of ", num_of_plates)
      plates_with_all_samps_per_exp_grp_to_print <- paste0(plates_with_all_samps_per_exp_grp, " out of ", num_of_plates)
      
      row_to_add <- c(batch_col, plates_with_no_samps_per_exp_grp_to_print, plates_with_all_samps_per_exp_grp_to_print)
      
      balloon_table[balloon_table_index,] <- row_to_add
    }
    return(balloon_table)
  })
  
  output$balloon_table <- renderTable({
    balloon_table_reac()
  })
  
  # take aways from distri of samples and missingness
  # Note after pareto about missinness in the whole data
  output$distri_text <- renderText({paste("1) An equal distribution of samples across batches, specifically the batch to correct on, is necessary for accurate batch correction.", "<br>", "2) Use the variation in sample distribution column from the above table to take a call on whether the data uploaded is suitable for batch correction.", "<br>", "3) Variation is calculated using the formula: (plate with max number of samples - plate with min number of samples/plate with max number of samples)*100, and we suggest having < 50% variation in sample distribution.", "<br>", "4) It is also recommended to have atleast 25 samples per plate and the \"Num_of_plates_passing_distri_cutoff\" column shows the number of plates passing this criteria.")})
  output$missing_text <- renderText({paste("1) Patterns in missingness distribution can make the data difficult to correct on.", "<br>", "2) Use the missingness plots and table above to ensure that batch correction can be performed accurately on your data.", "<br>", "3) Missingness should ideally be < 50% within each plate, the \"Num_of_plates_passing_missing_cutoff\" column shows how many plates passed this criteria.", "<br>", "4) Variation in missingness is calculated using the formula: (plate with max % of missingness - plate with min % of missingness/plate with max % of missingness).", "<br>", "The variation in missingness column tells you if there are plates with an uneven missingness distribution. This value should ideally be < 50%, which suggests that missingness is random and not particular to any plate!")})
  output$missing_note1 <- renderText({"As a note:"})
  output$missing_note2 <- renderText({paste("Count of missingness within your entire dataset (i.e., cells in the matrix with NAs) is ", sum(joined_df()$value), "and the percent of missing data is ", round((sum(joined_df()$value)/(nrow(sample_matrix())*(ncol(sample_matrix())-1)))*100, 2), "%. This percentage of missingness in overall data should also be < 50% for effective batch correction.")})
  output$balloon_text <- renderText({paste("1) It is recommended to have samples from each experimental group distributed across plates, without having plates with zero samples or all samples from one experimental group concentrated in one plate.", "<br>", "2) From the table above, it is ideal to have zero plates falling under either of the criteria.", "<br>", "3) Refer to the sample matrix plots to check for plates that have an unqueal sample distriution.")})
  output$take_away_text <- renderText({"If certain samples/batches has more missingness than others, filtering for fragments and samples can be applied in the next step. If you want to proceed with filtering and batch correction using the uploaded data, proceed to the Diagnosis+Filtering tab."})
  
  ############### Diagnosis+Filtering tab
  ########### Process samps_for_corr
  #combos <- reactive(NULL)
  observe ({
    #req(input$cols_of_int)
    anno_file_temp <- annotation_data()
    exp_grp <- exp_grp_val()
    unique_exp_grp <- unique(anno_file_temp[[exp_grp]])
    print(unique_exp_grp)
    cols_of_int_temp <- rv()
    print(cols_of_int_temp)
    combos <- tidyr::expand_grid(unique_exp_grp, cols_of_int_temp)
    combos <- combos %>%
      unite("combo", unique_exp_grp:cols_of_int_temp, sep= ":", 
            remove = FALSE)
    
    updateSelectInput(session, "samps_for_corr",
                      choices = combos$combo)
  })
  
  rv2 <- reactiveVal(NULL)
  observeEvent(input$samps_for_corr, {
    rv2(input$samps_for_corr)
  })
  
  samps_for_corr_val <- reactive({
    req(input$samps_for_corr)
    val <- as.character(rv2())
    val <- paste(as.character(val), collapse=", ")
    #val <- dQuote(val)
    val <- paste0("\"",val,"\"")
    val <- gsub(" ", "", val)
    return(val)
  })
  
  output_prefix_val <- reactive({ 
    val <- input$output_prefix
    val <- paste0("\"",val,"\"")
    return(val)
  })
  
  ########### Process iRT protein name
  iRT_prot_val <- reactive({ 
    val <- input$iRT_prot
    val <- paste0("\"",val,"\"")
    return(val)
  })
  
  ############### Calculate pvca and get var to correct on
  data_for_init_pvca <- reactive({
    req(input$submitId2)
    
    sample_matrix <- sample_matrix()
    protein <- protein_val()
    rownames(sample_matrix) <- sample_matrix[[protein]]
    
    sample_matrix <- sample_matrix[,-which(colnames(sample_matrix) %in% c(protein))] 
    sample_norm_long <- matrix_to_long(sample_matrix) 
    sample_norm_long.logged <- log_transform_df(sample_norm_long, log_base = 2, offset = 1)
    
    sample_matrix.logged <- long_to_matrix(sample_norm_long.logged) # convert back from long to matrix format
    sample_matrix.imputed_logged <- nafunctions(sample_matrix.logged, input$impute_method) # impute
    sample_long.logged <- matrix_to_long(sample_matrix.logged) # convert to long again, needed for batch correction
    
    return(sample_matrix.imputed_logged)
  })
  
  ################# reactive function that returns plot, data table and variable to correct on. 
  pvca_before_bc_reac <- reactive({
    data_for_init_pvca <- data_for_init_pvca()
    annotation_data <- annotation_data()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    pvca_before <- plot_PVCA(data_for_init_pvca,
                             annotation_data,
                             technical_factors = technical_factors,
                             biological_factors = biological_factors,
                             plot_title = "Before Batch Correction")
    plot <- pvca_before
    
    pvca_plot_data <- ggplot_build(plot)$plot$data
    var_to_correct_on <- toString(pvca_plot_data[1, "label"])
    
    pvca_plot <- annotate_figure(plot, top = text_grob("Before Batch Correction", color = "black", face = "bold", size = 15))
    
    to_return <- list(pvca_plot_data = pvca_plot_data, var_to_correct_on = var_to_correct_on, pvca_plot = pvca_plot)
    return(to_return)
  })
  
  output$results_msg <- renderText({
    to_return <- "PVCA is plotted here to find the batch/column to correct on. More input options will be visible on the left side panel once PVCA is plotted."
  })
  
  output$pvca_before_bc_init <- renderPlot({
    data <- pvca_before_bc_reac()
    pvca_plot <- data$pvca_plot
    print(pvca_plot)
  })
  
  pvca_plot_data_table <- reactive ({
    data <- pvca_before_bc_reac()
    data_table <- data$pvca_plot_data
    return(data_table)
  })
  
  var_to_correct_on <- reactive ({
    data <- pvca_before_bc_reac()
    var_to_correct_on <- data$var_to_correct_on
    return(var_to_correct_on)
  })
  
  output$pvca_table <- renderText({
    var_to_correct_on <- var_to_correct_on()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    exp_grp <- exp_grp_val()
    
    if (var_to_correct_on %in% check_cols_of_interest) {
      to_return <- paste0("NOTE: Batch-effect is most seen on ", var_to_correct_on, ", so batch correction will happen based on this column. If you want to change the column to correct on, choose the corresponsing variable in the drop-down menu in the side panel and press continue. If you are OK to correct on the chosen variable, directly press continue.")
      return(to_return) 
    }
    else if (var_to_correct_on == exp_grp) {
      to_return <- paste0("NOTE: Batch-effect is most seen on ", var_to_correct_on, ", so batch correction is not required on your data since most variation is seen in experimental group. If you still want to correct on one of the batches, choose the variable to correct on in the left side panel and press continue.")
      return(to_return)
    }
    else {
      to_return <- paste0("NOTE: Batch-effect is most seen on ", var_to_correct_on, ", so batch correction cannot be performed on your data since variation is seen from multiple factors in your data. If you still want to correct on one of the batches, choose the variable to correct on in the left side panel and press continue.")
      return(to_return)
    }
  })
  
  ########### Filtering input and filtering after pvca 
  ########### Process sample_threshold
  sample_threshold_val <- reactive({ 
    val <- input$sample_threshold
    return(val)
  })
  
  ########### Process expgrp_threshold
  expgrp_threshold_val <- reactive({ 
    val <- input$expgrp_threshold
    return(val)
  })
  
  ########### Process batch_threshold
  batch_threshold_val <- reactive({ 
    val <- input$batch_threshold
    return(val)
  })
  
  var_to_correct_on_final <- reactive ({
    #data <- pvca_plot_data()
    #var_to_correct_on <- data$var_to_correct_on
    var_to_correct_on_final <- input$var_to_correct_on
    return(var_to_correct_on_final)
  })
  
  observe({
    shinyjs::hide("output_prefix")
    shinyjs::hide("iRT_prot")
    shinyjs::hide("samps_for_corr")
    shinyjs::hide("var_to_correct_on")
    shinyjs::hide("batch_threshold")
    shinyjs::hide("expgrp_threshold")
    shinyjs::hide("sample_threshold")
    shinyjs::hide("continue")
    shinyjs::hide("results")
    if (!is.null(var_to_correct_on())) {
      shinyjs::show("output_prefix")
      shinyjs::show("iRT_prot")
      shinyjs::show("samps_for_corr")
      shinyjs::show("var_to_correct_on")
      shinyjs::show("batch_threshold")
      shinyjs::show("expgrp_threshold")
      shinyjs::show("sample_threshold")
      shinyjs::show("continue")
      shinyjs::disable("submitId2")
      updateSelectInput(session, "var_to_correct_on",
                        choices = unlist(strsplit(as.character(rv()), " ")),
                        selected = var_to_correct_on())
    }
  })
  
  observeEvent(input$submitId2, {
    # As soon as submit is hit, settings and impute method is disabled
    shinyjs::disable("data_type")
    shinyjs::disable("norm_file")
    shinyjs::disable("unnorm_file")
    shinyjs::disable("protein")
    shinyjs::disable("anno_file")
    shinyjs::disable("cols_of_int")
    shinyjs::disable("exp_grp")
    shinyjs::disable("sample_names")
    shinyjs::disable("impute_method")
  })
  
  observeEvent(input$continue, {
    # Show the plots and data before filt
    shinyjs::show("line_pvca")
    shinyjs::show("heading_filt")
    shinyjs::show("heading_perc_missing_plot")
    shinyjs::show("line_filt_plot1")
    shinyjs::show("heading_missing_expgrp")
    shinyjs::show("line_filt_plot2")
    shinyjs::show("heading_missing_batch")
    
    # hide other buttons and display new Go to Results button, also enable results tab
    shinyjs::hide("continue")
    shinyjs::hide("submitId2")
    shinyjs::show("results")
    shinyjs::enable(selector = '.navbar-nav a[data-value="Results"')
    
    # disable things in Diag+filt page
    shinyjs::disable("output_prefix")
    shinyjs::disable("iRT_prot")
    shinyjs::disable("samps_for_corr")
    shinyjs::disable("var_to_correct_on")
    shinyjs::disable("sample_names")
    shinyjs::disable("sample_threshold")
    shinyjs::disable("expgrp_threshold")
    shinyjs::disable("batch_threshold")
  })
  
  observeEvent(input$results, {
    updateNavbarPage(
      inputId = "maintab",
      selected = "Results"
    )
  })
  
  ############## Plots showing where the thresholds are and how much might get filtered out
  filt_plot1_reac <- reactive({
    req(input$continue)
    
    num_of_file <- length(rownames(sample_matrix()))
    if ((10/num_of_file)*40 > 10){
      FONTSIZE_sample_axis = 10
    }else{
      FONTSIZE_sample_axis = (10/num_of_file)*150
    }
    
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix(), annotation_data(), exp_grp, protein)
    filt_plot1 <- plot_per_sample(overlaps, FONTSIZE_sample_axis, sample_threshold_val(), color_list, exp_grp)
    return(filt_plot1)
  })
  
  output$filt_plot1 <- renderPlot({
    p <- filt_plot1_reac()
    print(p)
  })
  
  filt_plot2_reac <- reactive({
    req(input$continue)
    
    num_of_file <- length(rownames(sample_matrix()))
    if ((10/num_of_file)*40 > 10){
      FONTSIZE_sample_axis = 10
    }else{
      FONTSIZE_sample_axis = (10/num_of_file)*150
    }
    
    exp_grp <- exp_grp_val()
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix(), annotation_data(), exp_grp, protein)
    filt_plot2 <- plot_missgroup(overlaps, expgrp_threshold_val(), protein)
    return(filt_plot2)
  })
  
  output$filt_plot2 <- renderPlot({
    p <- filt_plot2_reac()
    print(p)
  })
  
  filt_plot3_reac <- reactive({
    req(input$continue)
    
    num_of_file <- length(rownames(sample_matrix()))
    if ((10/num_of_file)*40 > 10){
      FONTSIZE_sample_axis = 10
    }else{
      FONTSIZE_sample_axis = (10/num_of_file)*150
    }
    
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix(), annotation_data(), var_to_correct_on_final(), protein)
    filt_plot3 <- plot_missbatch(overlaps, batch_threshold_val(), protein)
    return(filt_plot3)
  })
  
  output$filt_plot3 <- renderPlot({
    p <- filt_plot3_reac()
    print(p)
  })
  
  ###### Filter samples with minimal and threshold filtering
  combo_output <- reactive({
    req(input$continue)
    
    sample_matrix <- sample_matrix()
    sample_matrix <- minimal_filtering(sample_matrix, annotation_data(), var_to_correct_on_final(), protein_val())
    sample_matrix <- filter_samples(sample_matrix, annotation_data(), sample_threshold_val(), expgrp_threshold_val(), var_to_correct_on_final(), batch_threshold_val(), exp_grp_val(), protein_val())
    
    sample_matrix_unnorm <- sample_matrix_unnorm()
    sample_matrix_unnorm <- minimal_filtering(sample_matrix_unnorm, annotation_data(), var_to_correct_on_final(),  protein_val())
    sample_matrix_unnorm <- filter_samples(sample_matrix_unnorm, annotation_data(), sample_threshold_val(), expgrp_threshold_val(), var_to_correct_on_final(), batch_threshold_val(),exp_grp_val(), protein_val())
    
    annotation_data <- as.data.frame(sample_matrix[[2]])
    sample_matrix <- as.data.frame(sample_matrix[[1]])
    sample_matrix_unnorm <- as.data.frame(sample_matrix_unnorm[[1]])
    
    combo <- list(sample_matrix = sample_matrix, sample_matrix_unnorm = sample_matrix_unnorm, annotation_data = annotation_data)
    
    return(combo)
  })
  
  sample_matrix_filt <- reactive({
    combo <- combo_output()
    sample_matrix_filt <- combo$sample_matrix
    return(sample_matrix_filt)
  })
  
  sample_matrix_unnorm_filt <- reactive({
    combo <- combo_output()
    sample_matrix_unnorm_filt <- combo$sample_matrix_unnorm
    return(sample_matrix_unnorm_filt)
  })
  
  annotation_data_filt <- reactive({
    combo <- combo_output()
    annotation_data_filt <- combo$annotation_data
    return(annotation_data_filt)
  })
  
  output$init_samps <- renderText({
    init_samps <- paste("Initial number of samples was: ", ncol(sample_matrix()))
    return(init_samps)
  })
  
  output$init_frags <- renderText({
    init_frags <- paste("Initial number of fragments was: ", nrow(sample_matrix()))
    return(init_frags)
  })
  
  output$filt_samps <- renderText({
    filt_samps <- paste("Number of samples remaining after filtering are: ", ncol(sample_matrix_filt()))
    return(filt_samps)
  })
  
  output$filt_frags <- renderText({
    filt_frags <- paste("Number of fragments remaining after filtering are: ", nrow(sample_matrix_filt()))
    return(filt_frags)
  })
  
  # For missingness in samples
  # Make table with sample_matrix and annotation combined
  joined_df2 <- reactive({
    req(input$cols_of_int)
    #req(input$unnorm_file, input$cols_of_int)
    table2 <- sample_matrix_filt() %>%
      select(everything()) %>%
      summarise_all(funs(sum(is.na(.)))) %>%
      rownames_to_column() %>%
      pivot_longer(cols=-rowname)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    
    joined_df2 <- merge(table2, annotation_data_filt(), by.x = "name", by.y = "FullRunName")
    joined_df2 <- joined_df2 %>% 
      select(check_cols_of_interest, exp_grp_val(), "name", "value")
    
    return(joined_df2)
  })
  
  output$filt_stats_note1 <- renderText({paste("Count of missingness within your entire dataset (i.e., cells in the matrix with NAs) was ", sum(joined_df()$value), "and the percent of missing data was ", round((sum(joined_df()$value)/(nrow(sample_matrix())*(ncol(sample_matrix())-1)))*100, 2), "%.")})
  output$filt_stats_note2 <- renderText({paste("Count of missingness within your entire dataset (i.e., cells in the matrix with NAs) after filtering is ", sum(joined_df2()$value), "and the percent of missing data is ", round((sum(joined_df2()$value)/(nrow(sample_matrix_filt())*(ncol(sample_matrix_filt())-1)))*100, 2), "%. This is the amout of data that will be imputed for calculating/plotting PCA and PVCA.")})
  
  ########### Generate plots with filtered results
  filtered_plot1_reac <- reactive({
    req(input$cols_of_int, input$exp_grp, input$sample_names)
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    technical_factors <- as.character(rv())
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors, technical_factors)
    print(selected_annotations)
    
    index = 0
    p_list <- list()
    for (i in selected_annotations) {
      index = index + 1
      print(i)
      print(index)
      #cat("### Samples by ", i, "\n\n", sep="")
      p_list[[index]] <- plot_bar(annotation_data_filt(), i, color_list)
      #print(p)
    }
    
    return(grid.arrange(grobs=p_list, ncol=3))
  })
  
  output$filtered_plot1 <- renderPlot({
    p <- filtered_plot1_reac()
    print(p)
  }, height = 400, width = 1000)
  
  # Balloon plot
  filt_balloon_plots_reac <- reactive({
    filt_balloon_plot_output_list <- lapply(1:length(as.character(rv())), function(i) {
      plotname <- paste("filt_balloonplot", i, sep="")
      plotOutput(plotname, height=500)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items to display properly.
    do.call(tagList, filt_balloon_plot_output_list)
  })
  
  num_plots <- 5
  for (i in 1:num_plots) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("filt_balloonplot", my_i, sep="")
      output[[plotname]] <- renderPlot({
        req(input$cols_of_int)
        
        index <- as.list(rv()[my_i])
        index_char <- toString(index)
        print(index_char)
        
        exp_grp <- exp_grp_val()
        balloon_plot <- plot_balloon(annotation_data_filt(), index_char, exp_grp)
        return(balloon_plot)
      })
    })
  }
  
  output$filtered_balloon_plots <- renderUI({
    p <- filt_balloon_plots_reac()
    print(p)
  })
  
  balloon_plot_reac_report2 <- reactive({
    req(input$cols_of_int)
    
    balloon_plots <- c()
    
    technical_factors <- unlist(strsplit(as.character(rv()), " "))
    
    for (i in 1:length(technical_factors)) {
      
      balloon_plots[[i]] <- local({
        i <- i
        index_char <- technical_factors[[i]]
        print(index_char)
        
        exp_grp <- exp_grp_val()
        balloon_plot <- plot_balloon(annotation_data_filt(), index_char, exp_grp)
        print(balloon_plot)
      })
    }
    return(balloon_plots)
  })
  
  # For missingness in samples
  # Pareto plot
  filtered_plot2_reac <- reactive({
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix_filt(), annotation_data_filt(), var_to_correct_on_final(), protein)
    plot <- plot_pareto(overlaps, protein)
    return(plot)
  })
  
  output$filtered_plot2 <- renderPlot({
    p <- filtered_plot2_reac()
    print(p)
  })
  
  # Heatmaps showing missingness
  # Insert the right number of plot output objects into the web page
  filtered_plots3_reac <- reactive({
    plot_output_list3 <- lapply(1:length(as.character(rv())), function(i) {
      plotname <- paste("missplot2", i, sep="")
      plotOutput(plotname)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items to display properly.
    do.call(tagList, plot_output_list3)
  })
  
  num_plots <- 5
  for (i in 1:num_plots) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("missplot2", my_i, sep="")
      output[[plotname]] <- renderPlot({
        req(input$cols_of_int)
        
        index <- as.list(rv()[my_i])
        index_char <- toString(index)
        print(index_char)
        
        num_of_file <- length(rownames(sample_matrix_filt()))
        if ((10/num_of_file)*40 > 10){
          FONTSIZE_sample_axis = 10
        }else{
          FONTSIZE_sample_axis = (10/num_of_file)*150
        }
        
        exp_grp <- exp_grp_val()
        sample_names <- sample_names_val()
        color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
        protein <- protein_val()
        overlaps <- call_missing_vals(sample_matrix_filt(), annotation_data_filt(), index_char, protein)
        miss_plot <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, index_char, protein)
        return(miss_plot)
      })
    })
  }
  
  output$filtered_plots3 <- renderUI({
    p <- filtered_plots3_reac()
    print(p)
  })
  
  missing_plots2_reac_report2 <- reactive({
    req(input$cols_of_int)
    
    missing_plots <- vector('list')
    
    technical_factors <- unlist(strsplit(as.character(rv()), " "))
    
    for (i in 1:length(technical_factors)) {
      
      missing_plots[[i]] <- local({
        i <- i
        index_char <- technical_factors[[i]]
        print(index_char)
        
        num_of_file <- length(rownames(sample_matrix_filt()))
        if ((10/num_of_file)*40 > 10){
          FONTSIZE_sample_axis = 10
        }else{
          FONTSIZE_sample_axis = (10/num_of_file)*150
        }
        
        exp_grp <- exp_grp_val()
        sample_names <- sample_names_val()
        color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
        protein <- protein_val()
        overlaps <- call_missing_vals(sample_matrix_filt(), annotation_data_filt(), index_char, protein)
        miss_plot <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, index_char, protein)
        #print(miss_plot)
      })
    }
    return(missing_plots)
  })
  
  filtered_plot4_reac <- reactive({
    req(input$cols_of_int)
    
    num_of_file <- length(rownames(sample_matrix_filt()))
    if ((10/num_of_file)*40 > 10){
      FONTSIZE_sample_axis = 10
    }else{
      FONTSIZE_sample_axis = (10/num_of_file)*150
    }
    
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix_filt(), annotation_data_filt(), exp_grp, protein)
    miss_plot2 <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, exp_grp, protein)
    return(miss_plot2)
  })
  
  output$filtered_plot4 <- renderPlot({
    p <- filtered_plot4_reac()
    print(p)
  })
  
  ######### getting data for visualizing normalization
  combo_output2 <- reactive({
      sample_matrix <- sample_matrix_filt()
      protein <- protein_val()
      rownames(sample_matrix) <- sample_matrix[[protein]] 
      sample_matrix <- sample_matrix[,-which(colnames(sample_matrix) %in% c(protein))] 
      sample_norm_long <- matrix_to_long(sample_matrix) 
      sample_norm_long.logged <- log_transform_df(sample_norm_long, log_base = 2, offset = 1)
      
      sample_matrix_unnorm <- sample_matrix_unnorm_filt()
      rownames(sample_matrix_unnorm) <- sample_matrix_unnorm[[protein]]
      sample_matrix_unnorm <- sample_matrix_unnorm[,-which(colnames(sample_matrix_unnorm) %in% c(protein))] 
      sample_unnorm_long <- matrix_to_long(sample_matrix_unnorm) 
      sample_unnorm_long.logged <- log_transform_df(sample_unnorm_long, log_base = 2, offset = 1)
      
      combo <- list(sample_norm_long.logged = sample_norm_long.logged, sample_unnorm_long.logged = sample_unnorm_long.logged)
      
      return(combo) 
  })
  
  sample_norm_long.logged <- reactive({
    combo <- combo_output2()
    sample_norm_long.logged <- combo$sample_norm_long.logged
    return(sample_norm_long.logged)
  })
  
  sample_unnorm_long.logged <- reactive({
    combo <- combo_output2()
    sample_unnorm_long.logged <- combo$sample_unnorm_long.logged
    return(sample_unnorm_long.logged)
  })
  
  unnorm_box_reac <- reactive({
    sample_unnorm_long.logged <- sample_unnorm_long.logged()
    annotation_data <- annotation_data_filt()
    exp_grp <- exp_grp_val()
    # Merge with annotation to get batch annotation
    sample_unnorm_long.logged[[exp_grp]] <- annotation_data[[exp_grp]][match(sample_unnorm_long.logged$FullRunName, annotation_data$FullRunName)]
    #sample_unnorm_long.logged$attribute_ExperimentalGroup <- annotation_data[["attribute_ExperimentalGroup"]][match(sample_unnorm_long.logged$FullRunName, annotation_data$FullRunName)]
    # order by batch and factor to plot in batch order
    sample_unnorm_long.logged <- sample_unnorm_long.logged[order(sample_unnorm_long.logged[[exp_grp]]),]
    sample_unnorm_long.logged$FullRunName = factor(sample_unnorm_long.logged$FullRunName, levels=unique(sample_unnorm_long.logged$FullRunName))
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    print(plot_boxnorm(sample_unnorm_long.logged, color_list, exp_grp))
  })
  
  output$unnorm_box <- renderPlot({
    p <- unnorm_box_reac()
    print(p)
  })
  
  norm_box_reac <- reactive({
    sample_norm_long.logged <- sample_norm_long.logged()
    annotation_data <- annotation_data_filt()
    exp_grp <- exp_grp_val()
    # Merge with annotation to get batch annotation
    sample_norm_long.logged[[exp_grp]] <- annotation_data[[exp_grp]][match(sample_norm_long.logged$FullRunName, annotation_data$FullRunName)]
    #sample_norm_long.logged$attribute_ExperimentalGroup <- annotation_data[["attribute_ExperimentalGroup"]][match(sample_norm_long.logged$FullRunName, annotation_data$FullRunName)]
    # order by batch and factor to plot in batch order
    sample_norm_long.logged <- sample_norm_long.logged[order(sample_norm_long.logged[[exp_grp]]),]
    sample_norm_long.logged$FullRunName = factor(sample_norm_long.logged$FullRunName, levels=unique(sample_norm_long.logged$FullRunName))
    
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    print(plot_boxnorm(sample_norm_long.logged, color_list, exp_grp))
  })
  
  output$norm_box <- renderPlot({
    p <- norm_box_reac()
  })
  
  #### Batch correction related
  before_correction_imputed_tables <- reactive({
    req(input$continue, input$impute_method)
    sample_norm_long.logged <- sample_norm_long.logged()
    
    sample_matrix.logged <- long_to_matrix(sample_norm_long.logged) # convert back from long to matrix format
    sample_matrix.imputed_logged <- nafunctions(sample_matrix.logged, input$impute_method) # impute
    sample_long.logged <- matrix_to_long(sample_matrix.logged) # convert to long again, needed for batch correction
    
    to_return <- list(sample_matrix.imputed_logged = sample_matrix.imputed_logged, sample_matrix.logged = sample_matrix.logged, sample_long.logged = sample_long.logged)
    
    return(to_return)
  })
  
  sample_matrix.imputed_logged <- reactive({
    from_return <- before_correction_imputed_tables()
    sample_matrix.imputed_logged <- from_return$sample_matrix.imputed_logged
    return(sample_matrix.imputed_logged)
  })
  
  sample_matrix.logged <- reactive({
    from_return <- before_correction_imputed_tables()
    sample_matrix.logged <- from_return$sample_matrix.logged
    return(sample_matrix.logged)
  })
  
  sample_long.logged <- reactive({
    from_return <- before_correction_imputed_tables()
    sample_long.logged <- from_return$sample_long.logged
    return(sample_long.logged)
  })
  
  pvca_before_bc_reac2 <- reactive({
    sample_matrix.imputed_logged <- sample_matrix.imputed_logged()
    annotation_data <- annotation_data()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    pvca_before <- plot_PVCA(sample_matrix.imputed_logged,
                             annotation_data,
                             technical_factors = technical_factors,
                             biological_factors = biological_factors,
                             plot_title = "Before Batch Correction")
    plot <- pvca_before
    
    pvca_plot_data_before_bc <- ggplot_build(plot)$plot$data
    plot_data_string <- toString(pvca_plot_data_before_bc)
    var_to_correct_on <- toString(pvca_plot_data_before_bc[1, "label"])
    
    pvca_dict_before_bc <- c()
    for (i in 1:nrow(pvca_plot_data_before_bc)) {
      #print(i)
      label = toString(pvca_plot_data_before_bc[i, "label"])
      weight = as.numeric(pvca_plot_data_before_bc[i, "weights"])
      pvca_dict_before_bc[label] <- weight
    }
    
    pvca_plot <- annotate_figure(plot, top = text_grob("Before Batch Correction", color = "black", face = "bold", size = 15))
    
    to_return <- list(pvca_dict_before_bc = pvca_dict_before_bc, pvca_plot = pvca_plot)
    return(to_return)
  })
  
  pvca_dict_before_bc <- reactive ({
    data <- pvca_before_bc_reac2()
    data_dict <- data$pvca_dict_before_bc
    return(data_dict)
  })
  
  hca_before_bc_reac <- reactive({
    req(input$continue)
    sample_matrix.imputed_logged.m <- data.matrix(sample_matrix.imputed_logged())
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors, technical_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    plot <- plotly_heatmap(sample_matrix.imputed_logged.m, df_with_annotations, color_list)
    return(plot)
  })
  
  pca_before_bc_reac <- reactive({
    req(input$continue)
    
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors, technical_factors)
    plot_for_pca <- c(biological_factors, technical_factors)
  
    pltList <- list()
    for( i in plot_for_pca ){
    # Create plot name.
      pltName <- i
      pltList[[ pltName ]] <- plot_PCA(sample_matrix.imputed_logged(), annotation_data_filt(), color_by = i,
                                     plot_title = i, color_scheme = color_list[[i]], width = 3, height = 3, units = c("in"))
    }
    plot <- grid.arrange(grobs=pltList, ncol=2,  top = text_grob("Before Batch Correction", color = "black", face = "bold", size = 18))
    return(plot)
  })
  
  frag_anno_for_iRT <- reactive({
    req(input$continue, input$iRT_prot)
    
    protein <- protein_val()
    
    # To make iRT plots, create own iRT mapping file
    sample_matrix <- sample_matrix_filt()
    iRT_prot_val <- input$iRT_prot
      
    fragment_annotation <- dplyr::filter(sample_matrix, grepl(iRT_prot_val, sample_matrix[[protein]]))
      
    fragment_annotation <- fragment_annotation %>% select(protein)
    colnames(fragment_annotation) <- c("peptide_group_label")
      
    fragment_annotation <- fragment_annotation %>%
      mutate(!! protein := iRT_prot_val(),
                  .before = "peptide_group_label")
      
    return(fragment_annotation) 
  })
  
  observe({
    if(input_provided(input$iRT_prot)){
      showTab("ResultsNav", "iRT mapping")
    }
    else{
      hideTab("ResultsNav", "iRT mapping")
    }
  })
  
  irt_before_bc_reac <- reactive({
    req(input$continue, input$iRT_prot)
    
    protein <- protein_val()
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    var_to_correct_on <- var_to_correct_on_final()
      
    psi_before <- plot_spike_in(sample_long.logged(), annotation_data_filt(),
                                peptide_annotation = frag_anno_for_iRT(),
                                protein_col = protein, spike_ins = iRT_prot_val(),
                                plot_title = 'Fragments - Before Batch Correction',
                                color_by_batch = TRUE, batch_col = var_to_correct_on, color_scheme = color_list[[var_to_correct_on]])
    psi_before <- psi_before + theme(legend.position="right")
    print(psi_before) 
  })
  
  observe({
    if(input_provided(input$samps_for_corr)){
      showTab("ResultsNav", "Correlation")
    }
    else{
      hideTab("ResultsNav", "Correlation")
    }
  })
  
  corr_before_bc_reac <- reactive({
    req(input$continue, input$samps_for_corr)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors, technical_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
          
    sample_for_correlation <-  c()
    column_for_correlation <-  c()
    correlation_samples <- unlist(strsplit(rv2(), ","))
    if (length(correlation_samples) > 0) {
      for( k in correlation_samples ) {
        split_samples_column <- unlist(strsplit(k, ":"))
        sample_for_correlation <- c(sample_for_correlation, split_samples_column[1])
        column_for_correlation <- c(column_for_correlation, split_samples_column[2])
      }
    }
        
    sample_element <- sample_for_correlation[1]
    col_element <- column_for_correlation[1]
          
    dr_replicate_filenames <- annotation_data_filt() %>%
      dplyr::filter(get(exp_grp) %in% sample_element) %>%
      dplyr::arrange(!!exp_grp, !!as.symbol(col_element)) %>%
      pull(!!sym('FullRunName'))
          
    replicate_filenames = as.character(dr_replicate_filenames)
          
    anno_for_corr <- annotation_data_filt() %>%
      dplyr::filter(FullRunName %in% replicate_filenames) %>%
      dplyr::select(FullRunName, col_element) %>%
      dplyr::group_by(across(all_of(col_element)))  %>%
      dplyr::arrange(across(all_of(col_element)))
          
    p1 = plotly_corr(sample_matrix.imputed_logged(), anno_for_corr, col_element, replicate_filenames, color_list)
    print(p1)
  })
  
  ## Plot same things in results tab
  output$pvca_before_bc_note <- renderText({
    text <- "Note: PVCA here might look different from the one seen in the diagnosis tab because this step uses the filtered data!"
    return(text)
  })
  
  output$pvca_before_bc <- renderPlot({
    req(input$continue)
    data <- pvca_before_bc_reac2()
    pvca_plot <- data$pvca_plot
    return(pvca_plot)
  })
  
  output$pca_before_bc <- renderPlot({
    req(input$continue)
    p <- pca_before_bc_reac()
    print(p)
  })
  
  #output$hca_before_bc <- renderPlotly({
  #  req(input$continue)
  #  p <- hca_before_bc_reac()
  #  print(p)
  #})
  
  output$irt_before_bc <- renderPlot({
    req(input$continue)
    p <- irt_before_bc_reac()
    print(p)
  })
  
  output$corr_before_bc <- renderPlotly({
    req(input$continue)
    p <- corr_before_bc_reac()
    print(p)
  })
  
  # Across batch-correction and plot
  across_correction_imputed_tables <- reactive({
    req(input$continue)
    
    # filtered, long, log transformed sample_matrix and filtered annotation data
    sample_long.logged <- sample_long.logged()
    annotation_data <- annotation_data_filt()
    
    # Order FullRunName in annotation with samples in input file
    annotation_data <- order_samples(sample_matrix_filt(), annotation_data)
    
    # Run across batch only
    sample_long.logged_bc_across <- correct_batch_effects_df(
      df_long = sample_long.logged,
      batch_col = var_to_correct_on_final(),
      sample_annotation = annotation_data,
      discrete_func = 'ComBat',
      continuous_func = NULL, 
      abs_threshold = 5,
      pct_threshold = 0.20
    )
    
    sample_matrix.logged_bc_across <- long_to_matrix(sample_long.logged_bc_across) # convert back from long to matrix format
    
    sample_matrix.imputed_logged_bc_across <- nafunctions(sample_matrix.logged_bc_across, input$impute_method) # impute
    
    # return sample matrix imputed, not imputed but batch-corrected, sample long form
    to_return <- list(sample_matrix.imputed_logged_bc_across = sample_matrix.imputed_logged_bc_across, sample_matrix.logged_bc_across = sample_matrix.logged_bc_across, sample_long.logged_bc_across = sample_long.logged_bc_across)
    return(to_return)
  })
  
  sample_matrix.imputed_logged_bc_across <- reactive({
    from_return <- across_correction_imputed_tables()
    sample_matrix.imputed_logged_bc_across <- from_return$sample_matrix.imputed_logged_bc_across
    return(sample_matrix.imputed_logged_bc_across)
  })
  
  sample_matrix.logged_bc_across <- reactive({
    from_return <- across_correction_imputed_tables()
    sample_matrix.logged_bc_across <- from_return$sample_matrix.logged_bc_across
    return(sample_matrix.logged_bc_across)
  })
  
  sample_long.logged_bc_across <- reactive({
    from_return <- across_correction_imputed_tables()
    sample_long.logged_bc_across <- from_return$sample_long.logged_bc_across
    return(sample_long.logged_bc_across)
  })
  
  #Two-way batch correction and plots
  after_correction_imputed_tables <- reactive({
    req(input$continue)
    
    # filtered, long, log transformed sample_matrix and filtered annotation data
    sample_long.logged <- sample_long.logged()
    annotation_data <- annotation_data_filt()
    
    # Order FullRunName in annotation with samples in input file
    annotation_data <- order_samples(sample_matrix_filt(), annotation_data)
    
    # batch correct and store
    # Run both within and across batches
    sample_long.logged_bc <- correct_batch_effects_df(
      df_long = sample_long.logged,
      batch_col = var_to_correct_on_final(),
      sample_annotation = annotation_data,
      discrete_func = 'ComBat',
      continuous_func = 'loess_regression', 
      abs_threshold = 5,
      pct_threshold = 0.20
    )
    
    sample_matrix.logged_bc <- long_to_matrix(sample_long.logged_bc) # convert back from long to matrix format
    
    sample_matrix.imputed_logged_bc <- nafunctions(sample_matrix.logged_bc, input$impute_method) # impute
    
    # return sample matrix imputed, not imputed but batch-corrected, sample long form
    to_return <- list(sample_matrix.imputed_logged_bc = sample_matrix.imputed_logged_bc, sample_matrix.logged_bc = sample_matrix.logged_bc, sample_long.logged_bc = sample_long.logged_bc)
  }) 
  
  sample_matrix.imputed_logged_bc <- reactive({
    from_return <- after_correction_imputed_tables()
    sample_matrix.imputed_logged_bc <- from_return$sample_matrix.imputed_logged_bc
    return(sample_matrix.imputed_logged_bc)
  })
  
  sample_matrix.logged_bc <- reactive({
    from_return <- after_correction_imputed_tables()
    sample_matrix.logged_bc <- from_return$sample_matrix.logged_bc
    return(sample_matrix.logged_bc)
  })
  
  sample_long.logged_bc <- reactive({
    from_return <- after_correction_imputed_tables()
    sample_long.logged_bc <- from_return$sample_long.logged_bc
    return(sample_long.logged_bc)
  })
  
  pvca_after_across_bc_reac <- reactive({
    req(input$continue)
    sample_matrix.imputed_logged_bc <- sample_matrix.imputed_logged_bc_across()
    annotation_data <- annotation_data_filt()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    pvca_after <- plot_PVCA(sample_matrix.imputed_logged_bc,
                            annotation_data,
                            technical_factors = technical_factors,
                            biological_factors = biological_factors,
                            plot_title = "After combat only batch-correction")
    plot <- pvca_after
    
    pvca_plot_data_after_across <- ggplot_build(plot)$plot$data
    plot_data_string <- toString(pvca_plot_data_after_across)
    var_to_correct_on <- toString(pvca_plot_data_after_across[1, "label"])
    
    pvca_dict_after_across <- c()
    for (i in 1:nrow(pvca_plot_data_after_across)) {
      #print(i)
      label = toString(pvca_plot_data_after_across[i, "label"])
      weight = as.numeric(pvca_plot_data_after_across[i, "weights"])
      #print(label)
      #print(weight)
      pvca_dict_after_across[label] <- weight
    }
    
    pvca_plot <- annotate_figure(plot, top = text_grob("After combat only batch-correction", color = "black", face = "bold", size = 15))
    
    to_return <- list(pvca_dict_after_across = pvca_dict_after_across, pvca_plot = pvca_plot)
    return(to_return)
  })
  
  pvca_dict_after_across <- reactive ({
    data <- pvca_after_across_bc_reac()
    data_dict <- data$pvca_dict_after_across
    return(data_dict)
  })
  
  output$pvca_after_across_bc <- renderPlot({
    data <- pvca_after_across_bc_reac()
    pvca_plot <- data$pvca_plot
    return(pvca_plot)
  })
  
  pvca_after_bc_reac <- reactive({
    req(input$continue)
    sample_matrix.imputed_logged_bc <- sample_matrix.imputed_logged_bc()
    annotation_data <- annotation_data_filt()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    pvca_after <- plot_PVCA(sample_matrix.imputed_logged_bc,
                            annotation_data,
                            technical_factors = technical_factors,
                            biological_factors = biological_factors,
                            plot_title = "After combat and loess batch-correction")
    plot <- pvca_after
    
    pvca_plot_data_after_bc <- ggplot_build(plot)$plot$data
    plot_data_string <- toString(pvca_plot_data_after_bc)
    var_to_correct_on <- toString(pvca_plot_data_after_bc[1, "label"])
    
    pvca_dict_after_bc <- c()
    for (i in 1:nrow(pvca_plot_data_after_bc)) {
      label = toString(pvca_plot_data_after_bc[i, "label"])
      weight = as.numeric(pvca_plot_data_after_bc[i, "weights"])
      pvca_dict_after_bc[label] <- weight
    }
    
    pvca_plot <- annotate_figure(plot, top = text_grob("After combat and loess batch-correction", color = "black", face = "bold", size = 15))
    
    to_return <- list(pvca_dict_after_bc = pvca_dict_after_bc, pvca_plot = pvca_plot)
    return(to_return)
  })
  
  pvca_dict_after_bc <- reactive ({
    data <- pvca_after_bc_reac()
    data_dict <- data$pvca_dict_after_bc
    return(data_dict)
  })
  
  output$pvca_after_bc <- renderPlot({
    data <- pvca_after_bc_reac()
    pvca_plot <- data$pvca_plot
    return(pvca_plot)
  })
  
  pca_after_across_bc_reac <- reactive({
    req(input$continue)
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors, technical_factors)
    plot_for_pca <- c(biological_factors,technical_factors)
    
    pltList <- list()
    for( i in plot_for_pca ){
      # Create plot name.
      pltName <- i
      pltList[[ pltName ]] <- plot_PCA(sample_matrix.imputed_logged_bc_across(), annotation_data_filt(), color_by = i,
                                       plot_title = i, color_scheme = color_list[[i]], width = 3, height = 3, units = c("in"))
    }
    plot <- grid.arrange(grobs=pltList, ncol=2,  top = text_grob("After combat only batch-correction", color = "black", face = "bold", size = 18))
    return(plot)
  })
  
  output$pca_after_across_bc <- renderPlot({
    p <- pca_after_across_bc_reac()
    print(p)
  })
  
  pca_after_bc_reac <- reactive({
    req(input$continue)
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors, technical_factors)
    plot_for_pca <- c(biological_factors,technical_factors)
    
    pltList <- list()
    for( i in plot_for_pca ){
      # Create plot name.
      pltName <- i
      pltList[[ pltName ]] <- plot_PCA(sample_matrix.imputed_logged_bc(), annotation_data_filt(), color_by = i,
                                       plot_title = i, color_scheme = color_list[[i]], width = 3, height = 3, units = c("in"))
    }
    plot <- grid.arrange(grobs=pltList, ncol=2,  top = text_grob("After combat and loess batch-correction", color = "black", face = "bold", size = 18))
    return(plot)
  })
  
  output$pca_after_bc <- renderPlot({
    p <- pca_after_bc_reac()
    print(p)
  })
  
  hca_after_across_bc_reac <- reactive({
    req(input$continue)
    sample_matrix.imputed_logged.m <- data.matrix(sample_matrix.imputed_logged_bc_across())
    
    exp_grp <- exp_grp_val()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors,technical_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    plot <- plotly_heatmap(sample_matrix.imputed_logged.m, df_with_annotations, color_list)
    return(plot)
  })
  
  #output$hca_after_across_bc <- renderPlotly({
  #  p <- hca_after_across_bc_reac()
  #  print(p)
  #})
  
  hca_after_bc_reac <- reactive({
    req(input$continue)
    sample_matrix.imputed_logged.m <- data.matrix(sample_matrix.imputed_logged_bc())
    
    exp_grp <- exp_grp_val()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors,technical_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    plot <- plotly_heatmap(sample_matrix.imputed_logged.m, df_with_annotations, color_list)
    return(plot)
  })
  
  #output$hca_after_bc <- renderPlotly({
  #  p <- hca_after_bc_reac()
  #  print(p)
  #})
  
  irt_after_across_bc_reac <- reactive({
    req(input$continue, input$iRT_prot)
    
    protein <- protein_val()
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    var_to_correct_on <- var_to_correct_on_final()
    
    psi_before <- plot_spike_in(sample_long.logged_bc_across(), annotation_data_filt(),
                                peptide_annotation = frag_anno_for_iRT(),
                                protein_col = protein, spike_ins = iRT_prot_val(),
                                plot_title = 'Fragments - After Combat Only Batch Correction',
                                color_by_batch = TRUE, batch_col = var_to_correct_on, color_scheme = color_list[[var_to_correct_on]])
    psi_before <- psi_before + theme(legend.position="right")
    print(psi_before)
  })
  
  output$irt_after_across_bc <- renderPlot({
    p <- irt_after_across_bc_reac()
    print(p)
  })
  
  irt_after_bc_reac <- reactive({
    req(input$continue, input$iRT_prot)
    
    protein <- protein_val()
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    var_to_correct_on <- var_to_correct_on_final()
    
    psi_before <- plot_spike_in(sample_long.logged_bc(), annotation_data_filt(),
                                peptide_annotation = frag_anno_for_iRT(),
                                protein_col = protein, spike_ins = iRT_prot_val(),
                                plot_title = 'Fragments - After Combat and Loess Batch Correction',
                                color_by_batch = TRUE, batch_col = var_to_correct_on, color_scheme = color_list[[var_to_correct_on]])
    psi_before <- psi_before + theme(legend.position="right")
    print(psi_before)
  })
  
  output$irt_after_bc <- renderPlot({
    p <- irt_after_bc_reac()
    print(p)
  })
  
  corr_after_across_bc_reac <- reactive({
    req(input$continue, input$samps_for_corr)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors,technical_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    sample_for_correlation <-  c()
    column_for_correlation <-  c()
    correlation_samples <- unlist(strsplit(rv2(), ","))
    if (length(correlation_samples) > 0) {
      for( k in correlation_samples ) {
        split_samples_column <- unlist(strsplit(k, ":"))
        sample_for_correlation <- c(sample_for_correlation, split_samples_column[1])
        column_for_correlation <- c(column_for_correlation, split_samples_column[2])
      }
    }
    
    sample_element <- sample_for_correlation[1]
    col_element <- column_for_correlation[1]
    
    dr_replicate_filenames <- annotation_data_filt() %>%
      dplyr::filter(get(exp_grp) %in% sample_element) %>%
      dplyr::arrange(!!exp_grp, !!as.symbol(col_element)) %>%
      pull(!!sym('FullRunName'))
    
    replicate_filenames = as.character(dr_replicate_filenames)
    
    anno_for_corr <- annotation_data_filt() %>%
      dplyr::filter(FullRunName %in% replicate_filenames) %>%
      dplyr::select(FullRunName, col_element) %>%
      dplyr::group_by(across(all_of(col_element)))  %>%
      dplyr::arrange(across(all_of(col_element)))
    
    p1 = plotly_corr(sample_matrix.imputed_logged_bc_across(), anno_for_corr, col_element, replicate_filenames, color_list)
    print(p1)
  })
  
  output$corr_after_across_bc <- renderPlotly({
    p <- corr_after_across_bc_reac()
    print(p)
  })

  corr_after_bc_reac <- reactive({
    req(input$continue, input$samps_for_corr)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest    
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors,technical_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    sample_for_correlation <-  c()
    column_for_correlation <-  c()
    correlation_samples <- unlist(strsplit(rv2(), ","))
    if (length(correlation_samples) > 0) {
      for( k in correlation_samples ) {
        split_samples_column <- unlist(strsplit(k, ":"))
        sample_for_correlation <- c(sample_for_correlation, split_samples_column[1])
        column_for_correlation <- c(column_for_correlation, split_samples_column[2])
      }
    }
    
    sample_element <- sample_for_correlation[1]
    col_element <- column_for_correlation[1]
    
    dr_replicate_filenames <- annotation_data_filt() %>%
      dplyr::filter(get(exp_grp) %in% sample_element) %>%
      dplyr::arrange(!!exp_grp, !!as.symbol(col_element)) %>%
      pull(!!sym('FullRunName'))
    
    replicate_filenames = as.character(dr_replicate_filenames)
    
    anno_for_corr <- annotation_data_filt() %>%
      dplyr::filter(FullRunName %in% replicate_filenames) %>%
      dplyr::select(FullRunName, col_element) %>%
      dplyr::group_by(across(all_of(col_element)))  %>%
      dplyr::arrange(across(all_of(col_element)))
    
    p1 = plotly_corr(sample_matrix.imputed_logged_bc(), anno_for_corr, col_element, replicate_filenames, color_list)
    print(p1)
  })
  
  output$corr_after_bc <- renderPlotly({
    p <- corr_after_bc_reac()
    print(p)
  })
  
  ########## Downloads section
  # First figure out if combat only or combat+loess was required
  bc_type_selection <- reactive({
    req(input$continue)
    
    # get all the dictionaries containing pvca data
    pvca_dict_before_bc <- pvca_dict_before_bc()
    pvca_dict_after_across <- pvca_dict_after_across()
    pvca_dict_after_bc <- pvca_dict_after_bc()
    
    # get selected annotations
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(biological_factors,technical_factors)
    
    pvca_after_across_score <- 0
    for (i in selected_annotations) {
      # make sure exp grp and cols of int are in pvca data
      if (i %in% names(pvca_dict_after_across)) {
        # for exp grp, make sure after correction is more than half of before correction. Add score if it is. 
        if (i == exp_grp) {
          if (pvca_dict_after_across[i] > (pvca_dict_before_bc[i]/2)) {
            pvca_after_across_score = pvca_after_across_score + 1
          }
        }
        # for other cols of int make sure the variation is lesser than before correction
        else if (pvca_dict_after_across[i] < pvca_dict_before_bc[i]) {
          pvca_after_across_score = pvca_after_across_score + 1
        }
      }
      # when a column of int is not in pvca data
      else {
        # if exp grp is not in pvca data, that is bad because that means variation in exp grp is completely removed so make score zero
        if (i == exp_grp) {
          pvca_after_across_score = 0
        }
        # if other cols of int are not in pvca data, that means variation in it has disappeared after bc, which is good, so turn score to +1. 
        else {
          pvca_after_across_score = pvca_after_across_score + 1
        }
      }
    }
    print("pvca_after_across_score")
    print(pvca_after_across_score)
    
    pvca_after_bc_score <- 0
    for (i in selected_annotations) {
      if (i %in% names(pvca_dict_after_bc)) {
        if (i == exp_grp) {
          if (pvca_dict_after_bc[i] > (pvca_dict_before_bc[i]/2)) {
            pvca_after_bc_score = pvca_after_bc_score + 1
          }
        }
        else if (pvca_dict_after_bc[i] < pvca_dict_before_bc[i]) {
          pvca_after_bc_score = pvca_after_bc_score + 1
        }
      }
      else {
        if (i == exp_grp) {
          pvca_after_bc_score = 0
        }
        else {
          pvca_after_bc_score = pvca_after_bc_score + 1
        }
      }
    }
    print("pvca_after_bc_score")
    print(pvca_after_bc_score)
    
    # Get final score
    BATCH_COL <- var_to_correct_on_final()
    final_score <- ""
    
    # If the scores are same
    if (pvca_after_across_score == pvca_after_bc_score) {
      # when both scores for after across and after complete are in the pvca data
      # get the one which has least variation in pvca data and assing that as the final one to use
      if (BATCH_COL %in% names(pvca_dict_after_across) && BATCH_COL %in% names(pvca_dict_after_bc)) {
        after_across_bc_val <- pvca_dict_after_across[BATCH_COL]
        after_bc_val <-  pvca_dict_after_bc[BATCH_COL]
        if (after_bc_val < after_across_bc_val) {
          final_score = "complete"
        } else {
          final_score = "across"
        }
        # if only after across variation exists, it means after bc variation is completely removed, so choose after bc
      } else if (BATCH_COL %in% names(pvca_dict_after_across)) {
        final_score = "complete"
        # if only after bc variation exists, it means after across variation is completely removed, so choose after across
      } else {
        final_score = "across"
      }
      # if the scores are different, assign the greater one as the final score
    } else {
      if (pvca_after_bc_score > pvca_after_across_score) {
        final_score = "complete"
      } else {
        final_score = "across"
      }
    }
    
    return(final_score)
  })
  
  output$bc_type_selection_note <- renderText ({
    final_score <- bc_type_selection()
    if (final_score == "complete") {
      return("Based on PVCA plots, we recommend you use combat+loess (complete) data for downstream analysis since variation in experimental group is maintained, and variation in other columns of interest have reduced with this type of batch correction. These are the first two files in this section below.")
    } else if (final_score == "across") {
      return("Based on PVCA plots, we recommend you use combat only data for downstream analysis since variation in experimental group is maintained, and variation in other columns of interest have reduced with this type of batch correction. These are the last two files in this section below.")
    }
  })
  
  output$filtnorm <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_filtered_normalized.csv", sep = "")
    },
    content = function(file) {
      write.csv(sample_matrix_filt(), file, row.names = FALSE)
    }
  )
  
  output$filtunnorm <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_filtered_unnormalized.csv", sep = "")
    },
    content = function(file) {
      write.csv(sample_matrix_unnorm_filt(), file, row.names = FALSE)
    }
  )
  
  output$filtanno <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_filtered_anno.csv", sep = "")
    },
    content = function(file) {
      write.csv(annotation_data_filt(), file, row.names = FALSE)
    }
  )
  
  output$before_preimpute <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_before_preimpute.csv", sep = "")
    },
    content = function(file) {
      results_matrix <- unlog_dm(sample_matrix.logged(), log_base = 2, offset = 1) # unlog
      results_matrix <- cbind(rownames(results_matrix), results_matrix) # unset rownames as protein names (see above)
      colnames(results_matrix)[1] <- "Protein" # Bring back protein column (see above)
      write.csv(results_matrix, file, row.names = FALSE)
    }
  )
  
  output$after_combat_preimpute <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_after_combat_preimpute.csv", sep = "")
    },
    content = function(file) {
      results_matrix <- unlog_dm(sample_matrix.logged_bc_across(), log_base = 2, offset = 1) # unlog
      results_matrix <- cbind(rownames(results_matrix), results_matrix) # unset rownames as protein names (see above)
      colnames(results_matrix)[1] <- "Protein" # Bring back protein column (see above)
      write.csv(results_matrix, file, row.names = FALSE)
    }
  )
  
  output$after_preimpute <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_after_preimpute.csv", sep = "")
    },
    content = function(file) {
      results_matrix <- unlog_dm(sample_matrix.logged_bc(), log_base = 2, offset = 1) # unlog
      results_matrix <- cbind(rownames(results_matrix), results_matrix) # unset rownames as protein names (see above)
      colnames(results_matrix)[1] <- "Protein" # Bring back protein column (see above)
      write.csv(results_matrix, file, row.names = FALSE)
    }
  )
  
  output$before_postimpute <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_before_postimpute.csv", sep = "")
    },
    content = function(file) {
      results_matrix <- unlog_dm(sample_matrix.imputed_logged(), log_base = 2, offset = 1) # unlog
      results_matrix <- cbind(rownames(results_matrix), results_matrix) # unset rownames as protein names (see above)
      colnames(results_matrix)[1] <- "Protein" # Bring back protein column (see above)
      write.csv(results_matrix, file, row.names = FALSE)
    }
  )
  
  output$after_combat_postimpute <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_after_combat_postimpute.csv", sep = "")
    },
    content = function(file) {
      results_matrix <- unlog_dm(sample_matrix.imputed_logged_bc_across(), log_base = 2, offset = 1) # unlog
      results_matrix <- cbind(rownames(results_matrix), results_matrix) # unset rownames as protein names (see above)
      colnames(results_matrix)[1] <- "Protein" # Bring back protein column (see above)
      write.csv(results_matrix, file, row.names = FALSE)
    }
  )
  
  output$after_postimpute <- downloadHandler(
    filename = function() {
      paste(input$output_prefix, "_after_postimpute.csv", sep = "")
    },
    content = function(file) {
      results_matrix <- unlog_dm(sample_matrix.imputed_logged_bc(), log_base = 2, offset = 1) # unlog
      results_matrix <- cbind(rownames(results_matrix), results_matrix) # unset rownames as protein names (see above)
      colnames(results_matrix)[1] <- "Protein" # Bring back protein column (see above)
      write.csv(results_matrix, file, row.names = FALSE)
    }
  )
  
  output$report <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      withProgress(message = 'Rendering, please wait!', {
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
        # one list() container object for all the parameters
        # all the objects have unique names (keys)
        if(input_provided(input$iRT_prot) && input_provided(input$samps_for_corr)) {
          params <- list(input_plot1_reac = input_plot1_reac(),
                         balloon_plot_reac_report = balloon_plot_reac_report(),
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac_report = missing_plots2_reac_report(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         balloon_plot_reac_report2 = balloon_plot_reac_report2(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         missing_plots2_reac_report2 = missing_plots2_reac_report2(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac(),
                         #hca_before_bc_reac = hca_before_bc_reac(),
                         #hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         #hca_after_bc_reac = hca_after_bc_reac(),
                         irt_before_bc_reac = irt_before_bc_reac(),
                         irt_after_across_bc_reac = irt_after_across_bc_reac(),
                         irt_after_bc_reac = irt_after_bc_reac(),
                         corr_before_bc_reac = corr_before_bc_reac(),
                         corr_after_across_bc_reac = corr_after_across_bc_reac(),
                         corr_after_bc_reac = corr_after_bc_reac())
        } else if(input_provided(input$iRT_prot)) {
          params <- list(input_plot1_reac = input_plot1_reac(),
                         balloon_plot_reac_report = balloon_plot_reac_report(),
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac_report = missing_plots2_reac_report(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         balloon_plot_reac_report2 = balloon_plot_reac_report2(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         missing_plots2_reac_report2 = missing_plots2_reac_report2(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac(),
                         #hca_before_bc_reac = hca_before_bc_reac(),
                         #hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         #hca_after_bc_reac = hca_after_bc_reac(),
                         irt_before_bc_reac = irt_before_bc_reac(),
                         irt_after_across_bc_reac = irt_after_across_bc_reac(),
                         irt_after_bc_reac = irt_after_bc_reac())
        } else if(input_provided(input$samps_for_corr)) {
          params <- list(input_plot1_reac = input_plot1_reac(),
                         balloon_plot_reac_report = balloon_plot_reac_report(),
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac_report = missing_plots2_reac_report(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         balloon_plot_reac_report2 = balloon_plot_reac_report2(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         missing_plots2_reac_report2 = missing_plots2_reac_report2(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac(),
                         #hca_before_bc_reac = hca_before_bc_reac(),
                         #hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         #hca_after_bc_reac = hca_after_bc_reac(),
                         corr_before_bc_reac = corr_before_bc_reac(),
                         corr_after_across_bc_reac = corr_after_across_bc_reac(),
                         corr_after_bc_reac = corr_after_bc_reac())
        } else {
          params <- list(input_plot1_reac = input_plot1_reac(),
                         balloon_plot_reac_report = balloon_plot_reac_report(),
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac_report = missing_plots2_reac_report(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         balloon_plot_reac_report2 = balloon_plot_reac_report2(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         missing_plots2_reac_report2 = missing_plots2_reac_report2(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac())
                         #hca_before_bc_reac = hca_before_bc_reac(),
                         #hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         #hca_after_bc_reac = hca_after_bc_reac()
        }
  
        rmarkdown::render(tempReport, output_file = file,
                    params = params,
                    envir = new.env(parent = globalenv()))
      })
  })
}

shinyApp(ui = ui, server = server)