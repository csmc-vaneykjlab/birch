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

########## Functions from bc script
#Replace hyphens in colnames with underscore
colClean <- function(x){ 
  colnames(x) <- gsub("-", "_", colnames(x))
  colnames(x) <- gsub("\\.", "_", colnames(x))
  return(x)  
} 

#Concatenate protein, peptide and fragment into a single column called "Protein" in input file
protein_concat <- function(df){
  df$Protein <- paste(df$ProteinName, df$PeptideSequence, df$FragmentIon, sep="-")
  df = df[ , !(names(df) %in% c('ProteinName','PeptideSequence','FragmentIon','RT'))]
  df <- df %>% dplyr::select(Protein, everything()) #makes the protein column the first column
  return(df)
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
  anno_cols <- c(exp_grp, sample_names, "order")
  #df_cols <- c("ProteinName", "PeptideSequence", "FragmentIon")
  df_cols <- c("Protein")
  
  suppressMessages({
    technical_factors <- check_cols_of_interest
    #biological_factors <- c('attribute_ExperimentalGroup')
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
    plot_for_pca <- c(technical_factors,biological_factors)
    #plot_for_pca <- c(technical_factors,biological_factors,'order')
    
    color_list <- sample_annotation_to_colors(
      annotation,
      factor_columns = c(technical_factors,biological_factors),
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

# function to make pareto plot
plot_pareto <- function(overlaps) {
  # Get percent missingness per fragment
  missval <- overlaps %>%
    mutate(Intensity = ifelse(Intensity > 0, 1, 0)) %>% 
    dplyr::group_by(Protein) %>%
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
plot_missval <- function(overlaps, FONTSIZE_sample_axis, color_list, batch_col) {
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
    group_by(Protein) %>%
    mutate(Intensity = ifelse(Intensity > 0, TRUE, FALSE)) %>% 
    dplyr::filter(any(Intensity == F)) %>%
    dplyr::arrange(Group) %>%
    dplyr::select(-Group) %>%
    dplyr::mutate(Intensity = if_else(Intensity, 1, 0)) %>%
    pivot_wider(names_from = SampleName, values_from = Intensity) %>%
    column_to_rownames("Protein") %>%
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
  
  ggplot(stat2, aes(x = SampleName, y = sum, fill = Group)) +
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
  fil_df = fil_df[(fil_df$Protein %in% frags_to_keep_c),]
  
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
  overlaps <- call_missing_vals(fil_df, annotation, batch_type, "Protein")
  
  #rename protein column to "Protein"
  #names(overlaps)[names(overlaps) == protein] <- "Protein"
  
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
    tabsetPanel(
      id="maintab",
      tabPanel(
        "Settings",
        sidebarPanel(
          fileInput("norm_file", "1. Choose normalized file", accept = ".txt"),
          fileInput("unnorm_file", "2. Choose unnormalized file", accept = ".txt"),
          selectInput("protein", "3. Choose the column containing the protein name", choices=c(), multiple=FALSE),
          fileInput("anno_file", "4. Choose annotation file", accept = ".txt"),
          selectInput("cols_of_int", "5. Choose the columns to correct for batch-effect", choices=c(), multiple=TRUE),
          selectInput("exp_grp", "6. Choose the column containing the experimental group", choices=c(), multiple=FALSE),
          selectInput("sample_names", "7. Choose the column containing the sample names (note that this column should match the headers in normalized data file)", choices=c(), multiple=FALSE),
          actionButton("nextId", "Next", class = "btn-warning"),
        width = 3),
        mainPanel(
          h3("Input parameters"),
          p("On uploading your files, a preview will be available for review. The normalized and unnormalized file must contain a column with protein names and the remaining columns should be intensities corresponding to samples as listed in the annotation file. The anotation file should contain columns to correct for batch-effect, a biological experimental group in which you want to retain variation, and a column with sample names that match the normalized/unnormalized file. Each of these columns can be specified in the menu on the left panel after uploading the data files."),
          
          h5(strong("Normalized file preview:")),
          dataTableOutput("norm_table"),
          hr(style = "border-top: 1px solid #000000;"),
          
          h5(strong("Unnormalized file preview:")),
          dataTableOutput("unnorm_table"),
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
                    hr(style = "border-top: 1px solid #000000;"),
                   
                    h5(strong("Importance of distribution of samples across batches:")),
                    htmlOutput("distri_text"),
                    hr(style = "border-top: 1px solid #000000;"),
                   
                    h5(strong("Importance of missing data across samples and batches:")),
                    htmlOutput("missing_text"),
                    linebreaks(2),
                   
                    textOutput("missing_note1"),
                    textOutput("missing_note2"),
                    hr(style = "border-top: 1px solid #000000;"),
                   
                    h5(strong("What's next:")),
                    textOutput("take_away_text"), 
                    linebreaks(1)),
        ),
      ),
      tabPanel(
        "Diagnosis",
        sidebarPanel(
          radioButtons("impute_method", "Choose the method to be used for imputation", methods),
          textInput("output_prefix", label = "Pick a prefix for output files", value = "test"),
          textInput("iRT_prot", label = "Enter the name of your iRT protein as listed in your data file"),
          selectInput("samps_for_corr", "Choose the \"experimental group:batch\" combo for correlation. Pro tip: Use replicates here.", choices=c(), multiple=TRUE, selected=NULL),
          selectInput("var_to_correct_on", "Choose the variable you would like to correct on:", choices=c()),
          sliderInput(inputId = "sample_threshold", 
                      label = "Sample Threshold", 
                      value = 0.7, min = 0, max = 1),
          sliderInput(inputId = "expgrp_threshold", 
                      label = "Experimental group Threshold", 
                      value = 0.5, min = 0, max = 1),
          sliderInput(inputId = "batch_threshold", 
                      label = "Batch Threshold", 
                      value = 0.7, min = 0, max = 1),
          actionButton("submitId2", "Submit"),
          actionButton("continue", "Continue"),
          width = 3),
        mainPanel(
          h3("Diagnosis for batch to correct on"),
          textOutput("results_msg"),
          h5(strong("PVCA before correction and filtering")),
          plotOutput("pvca_before_bc_init") %>% withSpinner(color="#0dc5c1"),
          linebreaks(2),
          textOutput("pvca_table"),
          hr(style = "border-top: 1px solid #000000;"),
          
          h3("Diagnosis for filtering"),
          h5(strong("Plot showing percent of missingness within each sample")),
          plotOutput("filt_plot1")  %>% withSpinner(color="#0dc5c1"),
          hr(style = "border-top: 1px solid #000000;"),
          
          h5(strong("Missingness percent by experimental group")),
          plotOutput("filt_plot2")  %>% withSpinner(color="#0dc5c1"),
          hr(style = "border-top: 1px solid #000000;"),
          
          h5(strong("Missingness percent by batch to correct on")),
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
          tabPanel("HCA",
                   tabsetPanel(
                     tabPanel("Before correction", plotlyOutput("hca_before_bc") %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat only correction", plotlyOutput("hca_after_across_bc") %>% withSpinner(color="#0dc5c1")), 
                     tabPanel("After combat and loess correction", plotlyOutput("hca_after_bc") %>% withSpinner(color="#0dc5c1"))
                   )
                  ),
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
                   downloadButton("after_combat_preimpute", "After combat batch-correction, pre-imputed file"),
                   linebreaks(2),
                   downloadButton("after_preimpute", "After complete batch-correction, pre-imputed file"),
                   linebreaks(2),
                   downloadButton("after_combat_postimpute", "After combat batch-correction, post-imputation file"),
                   linebreaks(2),
                   downloadButton("after_postimpute", "After complete batch-correction, post-imputation file"), 
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
  # To validate errors
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
  shinyjs::disable(selector = '.nav-tabs a[data-value="Initial analysis"')
  shinyjs::disable(selector = '.nav-tabs a[data-value="Diagnosis"')
  shinyjs::disable(selector = '.nav-tabs a[data-value="Results"')
  
  ########### Display norm file
  output$norm_table <- renderDataTable({
    req(input$norm_file)
    file <- input$norm_file
    
    if (is.null(file)){
      returnValue(data.frame(Description="Please upload normalized file"))
    } else {
      data<-read.delim(file$datapath, header = TRUE)
      #if(c("ProteinName", "PeptideSequence", "FragmentIon") %in% names(data)) {
      return(datatable(data, options = list(pageLength = 10)))
      #} else{
      #  shinyalert("Column Error","Please upload normalized data file that contains columns/column-names with ProteinName, PeptideSequence and FragmentIon.",type="error")
      #  returnValue(data.frame(Error="Please fix column error and reload the file"))
      #}
    }
  })
  
  norm_file_name <- reactive({ 
    file <- input$norm_file
    return(file$datapath)
  })
  
  ########### Display unnorm file
  output$unnorm_table <- renderDataTable({
    req(input$unnorm_file)
    file <- input$unnorm_file
    
    if (is.null(file)){
      returnValue(data.frame(Description="Please upload unnormalized file"))
    } else {
      data<-read.delim(file$datapath, header = TRUE)
      #if(c("ProteinName", "PeptideSequence", "FragmentIon") %in% names(data)) {
      return(datatable(data, options = list(pageLength = 10)))
      #} else{
      #  shinyalert("Column Error","Please upload unnormalized data file that contains columns/column-names with ProteinName, PeptideSequence and FragmentIon.",type="error")
      #  returnValue(data.frame(Error="Please fix column error and reload the file"))
      #}
    }
  })
  
  unnorm_file_name <- reactive({ 
    file <- input$unnorm_file
    return(file$datapath)
  })
  
  ############### Process norm file
  sample_matrix_init <- reactive({
    req(input$norm_file)
    norm_dataframe <- input$norm_file
    norm_dataframe <- read.delim(norm_dataframe$datapath, header = TRUE)
    return(norm_dataframe)
  })
  
  ############# Process unnorm file
  sample_matrix_unnorm_init <- reactive({
    req(input$unnorm_file)
    unnorm_dataframe <- input$unnorm_file
    unnorm_dataframe <- read.delim(unnorm_dataframe$datapath, header = TRUE)
    return(unnorm_dataframe)
  })
  
  ########### Display anno file
  output$anno_table <- renderDataTable({
    req(input$anno_file)
    file <- input$anno_file
    
    if (is.null(file)){
      returnValue(data.frame(Description="Please upload annotation file"))
    } else {
      data<-read.delim(file$datapath, header = TRUE)
      #if(c("attribute_ExperimentalGroup", "Level3", "order") %in% names(data)) {
      return(datatable(data, options = list(pageLength = 10)))
      #} else{
        #shinyalert("Column Error","Please upload annotation file that contains columns/column names with attribute_ExperimentalGroup, Level3, and order.",type="error")
        #returnValue(data.frame(Error="Please fix column error and reload the file"))
      #}
    }
  })
  
  anno_file_name <- reactive({ 
    file <- input$anno_file
    return(file$datapath)
  })
  
  ########### Process anno file
  annotation_data_init <- reactive({
    req(input$anno_file)
    anno_dataframe <- input$anno_file
    anno_dataframe <- read.delim(anno_dataframe$datapath, header = TRUE)
    return(anno_dataframe)
  })
  
  ########### Process cols_of_int
  observe({
    updateSelectInput(session, "cols_of_int",
                      choices = colnames(annotation_data_init()),
                      selected = colnames(annotation_data_init()[0]))
  })
  
  rv <- reactiveVal(NULL)
  observeEvent(input$cols_of_int, {
    rv(input$cols_of_int)
  })
  
  ########## Process exp_grp
  observe({
    updateSelectInput(session, "exp_grp",
                      choices = colnames(annotation_data_init()),
                      selected = colnames(annotation_data_init()[0]))
  })
  
  exp_grp_val <- reactiveVal(NULL)
  observeEvent(input$exp_grp, {
    exp_grp_val(input$exp_grp)
  })
  
  ########## Process sample_names
  observe({
    updateSelectInput(session, "sample_names",
                      choices = colnames(annotation_data_init()),
                      selected = colnames(annotation_data_init()[0]))
  })
  
  sample_names_val <- reactiveVal(NULL)
  observeEvent(input$sample_names, {
    sample_names_val(input$sample_names)
  })
  
  ########## Process protein
  observe({
    updateSelectInput(session, "protein",
                      choices = colnames(sample_matrix_init()),
                      selected = colnames(sample_matrix_init()[0]))
  })
  
  protein_val <- reactiveVal(NULL)
  observeEvent(input$protein, {
    protein_val(input$protein)
  })
  
  ########### Process files - part 2 - calling functions to clean and check the data
  annotation_data <- reactive({
    req(input$anno_file, input$sample_names)
    anno_dataframe <- input$anno_file
    anno_dataframe <- read.delim(anno_dataframe$datapath, header = TRUE)
    sample_names <- sample_names_val()
    anno_dataframe <- colClean2(anno_dataframe, sample_names)
    anno_dataframe <- rename_level3(anno_dataframe, sample_names)
    return(anno_dataframe)
  })
  
  sample_matrix <- reactive({
    req(input$norm_file)
    norm_dataframe <- input$norm_file
    norm_dataframe <- read.delim(norm_dataframe$datapath, header = TRUE)
    norm_dataframe <- colClean(norm_dataframe)
    norm_dataframe = norm_dataframe[ , !(names(norm_dataframe) %in% c('RT'))]
    #norm_dataframe <- protein_concat(norm_dataframe)
    return(norm_dataframe)
  })
  
  sample_matrix_unnorm <- reactive({
    req(input$unnorm_file)
    unnorm_dataframe <- input$unnorm_file
    unnorm_dataframe <- read.delim(unnorm_dataframe$datapath, header = TRUE)
    unnorm_dataframe <- colClean(unnorm_dataframe)
    unnorm_dataframe = unnorm_dataframe[ , !(names(unnorm_dataframe) %in% c('RT'))]
    #unnorm_dataframe <- protein_concat(unnorm_dataframe)
    return(unnorm_dataframe)
  })
  
  
  ###### Activate initial anal and filtering tab
  observeEvent(input$nextId, {
    updateTabsetPanel(
      inputId = "maintab",
      selected = "Initial analysis"
    )
    shinyjs::enable(selector = '.nav-tabs a[data-value="Initial analysis"')
    shinyjs::enable(selector = '.nav-tabs a[data-value="Diagnosis"')
  })
  
  ########### Generate initial analysis results
  input_plot1_reac <- reactive({
    req(input$cols_of_int, input$exp_grp, input$sample_names)
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data(),  as.character(rv()), exp_grp, sample_names)
    technical_factors <- as.character(rv())
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
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
    req(input$norm_file, input$cols_of_int)
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
  
  # Pareto plot
  pareto_plot_reac <- reactive({
    # hard coding batch to first batch in cols of int 
    index <- as.list(rv()[1])
    index_char <- toString(index)
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix(), annotation_data(), index_char, protein)
    plot <- plot_pareto(overlaps)
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
        miss_plot <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, index_char)
        return(miss_plot)
      })
    })
  }
  
  output$missing_plots2 <- renderUI({
    p <- missing_plots2_reac()
    print(p)
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
    miss_plot2 <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, exp_grp)
    
    return(miss_plot2)
  })
  
  output$missing_plot3 <- renderPlot({
    p <- missing_plot3_reac()
    print(p)
  })
  
  init_table_reac <- reactive({
    req(input$cols_of_int)
    
    technical_factors <- unlist(strsplit(as.character(rv()), " "))
    init_table <- data.frame(matrix(ncol = 12, nrow=length(technical_factors)))
    x <- c("Batch_column", "Total_num_of_plates", "Max_num_of_samples", "Min_num_of_samples", "Variation_in_sample_distribution", "Plates_passing_distri_cutoff", "Num_of_plates_passing_distri_cutoff", "Max_missing_by_plate", "Min_missing_by_plate", "Variation_in_missingness", "Plates_passing_missing_cutoff", "Num_of_plates_passing_missing_cutoff")
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
      select(-Max_missing_by_plate, -Min_missing_by_plate, -Max_num_of_samples, -Min_num_of_samples, -Plates_passing_missing_cutoff, -Plates_passing_distri_cutoff, -Total_num_of_plates) %>%
      mutate(Variation_in_missingness = round(Variation_in_missingness, 4)) %>%
      mutate(Variation_in_sample_distribution = round(Variation_in_sample_distribution, 4))

    init_table$Variation_in_missingness <- paste0(init_table$Variation_in_missingness, " %")
    init_table$Variation_in_sample_distribution <- paste0(init_table$Variation_in_sample_distribution, " %")
    
    return(init_table)
  })
  
  output$init_table <- renderTable({
    init_table_reac()
  })
  
  # take aways from distri of samples and missingness
  # Note after pareto about missinness in the whole data
  output$distri_text <- renderText({paste("1) An equal distribution of samples across batches, specifically the batch to correct on, is necessary for accurate batch correction.", "<br>", "2) Use the variation in sample distribution column from the above table to take a call on whether the data uploaded is suitable for batch correction.", "<br>", "3) Variation is calculated using the formula: (plate with max number of samples - plate with min number of samples/plate with max number of samples)*100, and we suggest having < 50% variation in sample distribution.", "<br>", "4) It is also recommended to have atleast 25 samples per plate and the \"Num_of_plates_passing_distri_cutoff\" column shows the number of plates passing this criteria.")})
  output$missing_text <- renderText({paste("1) Patterns in missingness distribution can make the data difficult to correct on.", "<br>", "2) Use the missingness plots and table above to ensure that batch correction can be performed accurately on your data.", "<br>", "3) Missingness should ideally be < 50% within each plate, the \"Num_of_plates_passing_missing_cutoff\" column shows how many plates passed this criteria.", "<br>", "4) Variation in missingness is calculated using the formula: (plate with max % of missingness - plate with min % of missingness/plate with max % of missingness).", "<br>", "The variation in missingness column tells you if there are plates with an uneven missingness distribution. This value should ideally be < 50%, which suggests that missingness is random and not particular to any plate!")})
  output$missing_note1 <- renderText({"As a note:"})
  output$missing_note2 <- renderText({paste("Count of missingness within your entire dataset (i.e., cells in the matrix with NAs) is ", sum(joined_df()$value), "and the percent of missing data is ", round((sum(joined_df()$value)/(nrow(sample_matrix())*(ncol(sample_matrix())-1)))*100, 2), "%. This percentage of missingness in overall data should also be < 50% for effective batch correction.")})
  output$take_away_text <- renderText({"If certain samples/batches has more missingness than others, filtering for fragments and samples can be applied in the next step. If you want to proceed with filtering and batch correction using the uploaded data, proceed to the Diagnosis tab."})
  
  ############### Diagnosis tab
  ########### Process samps_for_corr
  #combos <- reactive(NULL)
  observe ({
    #req(input$cols_of_int)
    anno_file_temp <- annotation_data()
    unique_exp_grp <- unique(anno_file_temp$attribute_ExperimentalGroup)
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
    rownames(sample_matrix) <- sample_matrix$Protein 
    
    sample_matrix <- sample_matrix[,-which(colnames(sample_matrix) %in% c("Protein", "nPeptide", "nFragment", "RT"))] 
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
    shinyjs::hide("samps_for_corr")
    shinyjs::hide("var_to_correct_on")
    shinyjs::hide("batch_threshold")
    shinyjs::hide("expgrp_threshold")
    shinyjs::hide("sample_threshold")
    shinyjs::hide("continue")
    if (!is.null(var_to_correct_on())) {
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
    req(input$norm_file, input$cols_of_int)
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
    selected_annotations <- c(technical_factors,biological_factors)
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
  
  # For missingness in samples
  # Pareto plot
  filtered_plot2_reac <- reactive({
    protein <- protein_val()
    overlaps <- call_missing_vals(sample_matrix_filt(), annotation_data_filt(), var_to_correct_on_final(), protein)
    plot <- plot_pareto(overlaps)
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
        miss_plot <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, index_char)
        return(miss_plot)
      })
    })
  }
  
  output$filtered_plots3 <- renderUI({
    p <- filtered_plots3_reac()
    print(p)
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
    miss_plot2 <- plot_missval(overlaps, FONTSIZE_sample_axis, color_list, exp_grp)
    return(miss_plot2)
  })
  
  output$filtered_plot4 <- renderPlot({
    p <- filtered_plot4_reac()
    print(p)
  })
  
  ######### getting data for visualizing normalization
  combo_output2 <- reactive({
    sample_matrix <- sample_matrix_filt()
    rownames(sample_matrix) <- sample_matrix$Protein 
    sample_matrix <- sample_matrix[,-which(colnames(sample_matrix) %in% c("Protein", "nPeptide", "nFragment", "RT"))] 
    sample_norm_long <- matrix_to_long(sample_matrix) 
    sample_norm_long.logged <- log_transform_df(sample_norm_long, log_base = 2, offset = 1)
    
    sample_matrix_unnorm <- sample_matrix_unnorm_filt()
    rownames(sample_matrix_unnorm) <- sample_matrix_unnorm$Protein 
    sample_matrix_unnorm <- sample_matrix_unnorm[,-which(colnames(sample_matrix_unnorm) %in% c("Protein", "nPeptide", "nFragment", "RT"))] 
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
  
  ####### Activate Results tab
  observeEvent(input$continue, {
    shinyjs::enable(selector = '.nav-tabs a[data-value="Results"')
  })
  
  #observeEvent(input$continue, {
  #  updateTabsetPanel(
  #    inputId = "maintab",
  #    selected = "Results"
  #  )
  #})
  
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
    
    #pvca_plot_data <- ggplot_build(plot)$plot$data
    #plot_data_string <- toString(pvca_plot_data)
    #var_to_correct_on <- toString(pvca_plot_data[1, "label"])
    
    plot <- annotate_figure(plot, top = text_grob("Before Batch Correction", color = "black", face = "bold", size = 15))
    return(plot)
  })
  
  hca_before_bc_reac <- reactive({
    req(input$continue)
    sample_matrix.imputed_logged.m <- data.matrix(sample_matrix.imputed_logged())
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
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
    selected_annotations <- c(technical_factors,biological_factors)
    plot_for_pca <- c(technical_factors,biological_factors,'order')
  
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
    
    # To make iRT plots, create own iRT mapping file
    sample_matrix <- sample_matrix_filt()
    iRT_prot_val <- input$iRT_prot
      
    fragment_annotation <- dplyr::filter(sample_matrix, grepl(iRT_prot_val, sample_matrix$Protein))
      
    fragment_annotation <- fragment_annotation["Protein"]
    colnames(fragment_annotation) <- c("peptide_group_label")
      
    fragment_annotation <- fragment_annotation %>%
      add_column(Protein = iRT_prot_val(),
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
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)

    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    var_to_correct_on <- var_to_correct_on()
      
    psi_before <- plot_spike_in(sample_long.logged(), annotation_data_filt(),
                                peptide_annotation = frag_anno_for_iRT(),
                                protein_col = 'Protein', spike_ins = iRT_prot_val(),
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
    selected_annotations <- c(technical_factors,biological_factors)
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
      dplyr::filter(attribute_ExperimentalGroup %in% sample_element) %>%
      dplyr::arrange(attribute_ExperimentalGroup, !!as.symbol(col_element)) %>%
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
    p <- pvca_before_bc_reac2()
    print(p)
  })
  
  output$pca_before_bc <- renderPlot({
    req(input$continue)
    p <- pca_before_bc_reac()
    print(p)
  })
  
  output$hca_before_bc <- renderPlotly({
    req(input$continue)
    p <- hca_before_bc_reac()
    print(p)
  })
  
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
    plot <- annotate_figure(plot, top = text_grob("After combat only batch-correction", color = "black", face = "bold", size = 15))
    return(plot)
  })
  
  output$pvca_after_across_bc <- renderPlot({
    p <- pvca_after_across_bc_reac()
    print(p)
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
    plot <- annotate_figure(plot, top = text_grob("After combat and loess batch-correction", color = "black", face = "bold", size = 15))
    return(plot)
  })
  
  output$pvca_after_bc <- renderPlot({
    p <- pvca_after_bc_reac()
    print(p)
  })
  
  pca_after_across_bc_reac <- reactive({
    req(input$continue)
    exp_grp <- exp_grp_val()
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
    plot_for_pca <- c(technical_factors,biological_factors,'order')
    
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
    selected_annotations <- c(technical_factors,biological_factors)
    plot_for_pca <- c(technical_factors,biological_factors,'order')
    
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
    selected_annotations <- c(technical_factors,biological_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    plot <- plotly_heatmap(sample_matrix.imputed_logged.m, df_with_annotations, color_list)
    return(plot)
  })
  
  output$hca_after_across_bc <- renderPlotly({
    p <- hca_after_across_bc_reac()
    print(p)
  })
  
  hca_after_bc_reac <- reactive({
    req(input$continue)
    sample_matrix.imputed_logged.m <- data.matrix(sample_matrix.imputed_logged_bc())
    
    exp_grp <- exp_grp_val()
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    plot <- plotly_heatmap(sample_matrix.imputed_logged.m, df_with_annotations, color_list)
    return(plot)
  })
  
  output$hca_after_bc <- renderPlotly({
    p <- hca_after_bc_reac()
    print(p)
  })
  
  irt_after_across_bc_reac <- reactive({
    req(input$continue, input$iRT_prot)
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    var_to_correct_on <- var_to_correct_on_final()
    
    psi_before <- plot_spike_in(sample_long.logged_bc_across(), annotation_data_filt(),
                                peptide_annotation = frag_anno_for_iRT(),
                                protein_col = 'Protein', spike_ins = iRT_prot_val(),
                                plot_title = 'Fragments - After Combat Batch Correction',
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
    
    check_cols_of_interest <- unlist(strsplit(as.character(rv()), " "))
    technical_factors <- check_cols_of_interest
    exp_grp <- exp_grp_val()
    biological_factors <- c(exp_grp)
    selected_annotations <- c(technical_factors,biological_factors)
    df_with_annotations <- hca_anno_df(annotation_data_filt(), selected_annotations)
    sample_names <- sample_names_val()
    color_list <- create_color_list(annotation_data_filt(),  as.character(rv()), exp_grp, sample_names)
    
    var_to_correct_on <- var_to_correct_on_final()
    
    psi_before <- plot_spike_in(sample_long.logged_bc(), annotation_data_filt(),
                                peptide_annotation = frag_anno_for_iRT(),
                                protein_col = 'Protein', spike_ins = iRT_prot_val(),
                                plot_title = 'Fragments - After Batch Correction',
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
    selected_annotations <- c(technical_factors,biological_factors)
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
      dplyr::filter(attribute_ExperimentalGroup %in% sample_element) %>%
      dplyr::arrange(attribute_ExperimentalGroup, !!as.symbol(col_element)) %>%
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
    selected_annotations <- c(technical_factors,biological_factors)
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
      dplyr::filter(attribute_ExperimentalGroup %in% sample_element) %>%
      dplyr::arrange(attribute_ExperimentalGroup, !!as.symbol(col_element)) %>%
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
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac = missing_plots2_reac(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         filtered_plots3_reac = filtered_plots3_reac(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac(),
                         hca_before_bc_reac = hca_before_bc_reac(),
                         hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         hca_after_bc_reac = hca_after_bc_reac(),
                         irt_before_bc_reac = irt_before_bc_reac(),
                         irt_after_across_bc_reac = irt_after_across_bc_reac(),
                         irt_after_bc_reac = irt_after_bc_reac(),
                         corr_before_bc_reac = corr_before_bc_reac(),
                         corr_after_across_bc_reac = corr_after_across_bc_reac(),
                         corr_after_bc_reac = corr_after_bc_reac())
        } else if(input_provided(input$iRT_prot)) {
          params <- list(input_plot1_reac = input_plot1_reac(),
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac = missing_plots2_reac(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         filtered_plots3_reac = filtered_plots3_reac(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac(),
                         hca_before_bc_reac = hca_before_bc_reac(),
                         hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         hca_after_bc_reac = hca_after_bc_reac(),
                         irt_before_bc_reac = irt_before_bc_reac(),
                         irt_after_across_bc_reac = irt_after_across_bc_reac(),
                         irt_after_bc_reac = irt_after_bc_reac())
        } else if(input_provided(input$samps_for_corr)) {
          params <- list(input_plot1_reac = input_plot1_reac(),
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac = missing_plots2_reac(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         filtered_plots3_reac = filtered_plots3_reac(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac(),
                         hca_before_bc_reac = hca_before_bc_reac(),
                         hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         hca_after_bc_reac = hca_after_bc_reac(),
                         corr_before_bc_reac = corr_before_bc_reac(),
                         corr_after_across_bc_reac = corr_after_across_bc_reac(),
                         corr_after_bc_reac = corr_after_bc_reac())
        } else {
          params <- list(input_plot1_reac = input_plot1_reac(),
                         pareto_plot_reac = pareto_plot_reac(),
                         missing_plots2_reac = missing_plots2_reac(),
                         missing_plot3_reac = missing_plot3_reac(),
                         init_table_reac = init_table_reac(),
                         filt_plot1_reac = filt_plot1_reac(),
                         filt_plot2_reac = filt_plot2_reac(),
                         filt_plot3_reac = filt_plot3_reac(),
                         filtered_plot1_reac = filtered_plot1_reac(),
                         filtered_plot2_reac = filtered_plot2_reac(),
                         filtered_plots3_reac = filtered_plots3_reac(),
                         filtered_plot4_reac  = filtered_plot4_reac(),
                         unnorm_box_reac = unnorm_box_reac(),
                         norm_box_reac = norm_box_reac(),
                         pvca_before_bc_reac = pvca_before_bc_reac2(),
                         pvca_after_across_bc_reac = pvca_after_across_bc_reac(),
                         pvca_after_bc_reac = pvca_after_bc_reac(),
                         pca_before_bc_reac = pca_before_bc_reac(),
                         pca_after_across_bc_reac = pca_after_across_bc_reac(),
                         pca_after_bc_reac = pca_after_bc_reac(),
                         hca_before_bc_reac = hca_before_bc_reac(),
                         hca_after_across_bc_reac = hca_after_across_bc_reac(),
                         hca_after_bc_reac = hca_after_bc_reac())
        }
  
        rmarkdown::render(tempReport, output_file = file,
                    params = params,
                    envir = new.env(parent = globalenv()))
      })
  })
}

shinyApp(ui = ui, server = server)