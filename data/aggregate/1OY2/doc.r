# library(ggplot2)
# library(RColorBrewer)
# library(ggpubr)
# library(plotly)


ggplot_ethan_theme <- function(plot){
  
  plot + theme_minimal() + 
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      legend.text=element_text(size=10),
      legend.title=element_text(size=15),
      legend.key.height = unit(2.5, "cm"),
      axis.text = element_text(size = 14, face='bold')
    )
  
}

add_combined_metric <- function(agg.df){
  
  min_ts <- min(agg.df$total_score)
  max_ts <- max(agg.df$total_score)
  max_idx <- max(agg.df$interface_delta_X)
  min_idx <- min(agg.df$interface_delta_X)
  
  combined_list <- list()
  
  for (i in 1:nrow(agg.df)){
    ts <- agg.df[i, ]$total_score
    idx <- agg.df[i, ]$interface_delta_X
    cm <- 2 * ((ts - min_ts) / (max_ts - min_ts)) - 1 + 2*((idx - min_idx) / (max_idx - min_idx)) - 1
    combined_list[[i]] <- cm
  }
  agg.df$combined_metric <- unlist(combined_list)
  agg.df
  
}

energy_hole_plot_total_score <- function(agg.df, top_n=0.15){
  agg.df.ts <- agg.df[order(agg.df$total_score), ]
  agg.df.ts.top <- agg.df.ts[1:round(nrow(agg.df.ts) * top_n), ]
  p <- ggplot(agg.df.ts.top, aes(x=distance_to_pocket, y=total_score, color=interface_delta_X)) +
    geom_point(alpha=0.7, fill='black') +
    labs(x='Distance to PTB (Å)', y='Total Score', color='Interface Delta') +
    scale_colour_viridis_c()
  
  p <- ggplot_ethan_theme(p)
  p
}

energy_hole_plot_id <- function(agg.df, top_n=0.15){
  agg.df.ts <- agg.df[order(agg.df$interface_delta_X), ]
  agg.df.ts.top <- agg.df.ts[1:round(nrow(agg.df.ts) * top_n), ]
  p <- ggplot(agg.df.ts.top, aes(text=pdb_path, x=distance_to_pocket, y=interface_delta_X, color=total_score)) +
    geom_point(alpha=0.7, fill='black') +
    labs(x='Distance to PTB (Å)', y='Interface Delta', color='Total Score') +
    scale_colour_viridis_c()
  
  p <- ggplot_ethan_theme(p)
  p
}

write_csv <- function(df, file){
  
  write.csv(df, file, row.names = FALSE, quote = FALSE)
  
}

cut_and_subset <- function(df, file, select=c()){
  df.cut <- df[1:10, ]
  df.cut.subset <- subset(df.cut, select=select)
  df.cut.subset$pdb_path <- unlist(lapply(df.cut.subset$pdb_path, basename))
  write_csv(df.cut.subset, file)
}

cut_and_subset_move_to_target<- function(df, target_dir){
  dir.create(target_dir)
  df.cut <- df[1:10, ]
  df.cut.subset <- subset(df.cut, select='pdb_path')
  cmd_list <- list()
  i = 0
  for (pdb in df.cut.subset){
    i <- 1 + i
    # pdb <- gsub('farm_cluster/', 'farm_cluster/jobs/', pdb)
    cmd <- sprintf('cp %s %s', pdb, paste(target_dir, basename(pdb), sep='/'))
    cmd_list[[i]] <- cmd
  }
  print(cmd_list)
  lapply(cmd_list[[1]], system)
  cmd_list
}

add_combined_metric <- function(agg.df){
  
  min_ts <- min(agg.df$total_score)
  max_ts <- max(agg.df$total_score)
  max_idx <- max(agg.df$interface_delta_X)
  min_idx <- min(agg.df$interface_delta_X)
  
  combined_list <- list()
  
  for (i in 1:nrow(agg.df)){
    ts <- agg.df[i, ]$total_score
    idx <- agg.df[i, ]$interface_delta_X
    cm <- 2 * ((ts - min_ts) / (max_ts - min_ts)) - 1 + 2*((idx - min_idx) / (max_idx - min_idx)) - 1
    combined_list[[i]] <- cm
  }
  agg.df$combined_metric <- unlist(combined_list)
  agg.df
  
}

remote_pdb_download_top_10 <- function(df, temp_remote_dir){
  
  top_to_paths <- df[1:10, ]$pdb_path
  copy_command = "cp"
  for (p in top_to_paths){
    copy_command <- paste(copy_command, p)
  }
  copy_command <- paste(copy_command, temp_remote_dir)
  print(copy_command)
}

#agg.path <- "/home/ethan/share/research/Davis/rotations/Gino/random_docking_expanded_sampling/process_agg_results/high_aff_ligands_processed"
#agg.path <- "/home/ethan/share/backups/farm_cluster/jobs/trka_to_10y2_agg/trka_processed.rds"
#agg.path <- "/home/ethan/share/backups/farm_cluster/jobs/RTX73145433_agg/RTX73145433_agg.rds"
#agg.path <- '/home/ethan/share/backups/farm_cluster/jobs/trka_to_1shc/trka_1shc_processed.rds'
#agg.path <- '/home/ethan/Documents/github/docking_results/data/Trka/best_poses_tables/npey_agg.rds'
#agg.path <- '/home/ethan/Documents/github/docking_results/data/Trka/trka_agg.rds'
agg.path <- '/home/ethollem/jobs/random_docking_trka/aggregate/1OY2/1OY2_ag.rds'
agg.df <- readRDS(agg.path)
#agg.df <- as.data.frame(read.delim('/home/ethan/Documents/github/docking_results/data/Trka/best_poses_tables/npey_agg.tsv'))
agg.df <- add_combined_metric(agg.df)
agg.df.ts <- agg.df[order(agg.df$total_score), ]
agg.df.id <- agg.df[order(agg.df$interface_delta_X), ]
agg.df.c <- agg.df[order(agg.df$combined_metric), ]




#id.plt <- energy_hole_plot_id(agg.df)
#ts.plt <- energy_hole_plot_total_score(agg.df)
#g <- ggarrange(id.plt, ts.plt, nrow=1, ncol=2, labels=c('a', 'b'))


#remote_pdb_download_top_10(agg.df.c, '/home/ethollem/jobs/random_docking_trka/top_poses/combined_metric')
#remote_pdb_download_top_10(agg.df.id, '/home/ethollem/jobs/random_docking_trka/top_poses/total_score')
#remote_pdb_download_top_10(agg.df.ts, '/home/ethollem/jobs/random_docking_trka/top_poses/interface_delta_x')

pdb.dir <- '/home/ethollem/jobs/random_docking_trka/aggregate/1OY2/top_poses_pdb'
dir.create(pdb.dir)
pdb.dir.id <- file.path(pdb.dir, 'interface_delta_x')
dir.create(pdb.dir.id)
pdb.dir.ts <- file.path(pdb.dir, 'total_score')
dir.create(pdb.dir.ts)
pdb.dir.c <- file.path(pdb.dir, 'combined_metric')
dir.create(pdb.dir.c)

cut_and_subset_move_to_target(agg.df.id, pdb.dir.id)
cut_and_subset_move_to_target(agg.df.ts, pdb.dir.ts)
cut_and_subset_move_to_target(agg.df.c, pdb.dir.c)

tables_dir <- '/home/ethollem/jobs/random_docking_trka/aggregate/1OY2/tables'
dir.create(tables_dir)
id.table <- file.path(tables_dir, 'trka_1OY2_idx.csv')
ts.table <- file.path(tables_dir, 'trka_1OY2_ts.csv')
c.table <- file.path(tables_dir, 'trka_1OY2_cm.csv')

cut_and_subset(agg.df.id, id.table, c('pdb_path', 'interface_delta_X'))
cut_and_subset(agg.df.ts, ts.table, c('pdb_path', 'total_score'))
cut_and_subset(agg.df.c, c.table, c('pdb_path', 'combined_metric'))

# 
# 



#print(rbind(agg.df.id[1:4, ]$pdb_path, print(agg.df.id[1:4, ]$interface_delta_X)))
#print(rbind(agg.df.ts[1:4, ]$pdb_path, print(agg.df.ts[1:4, ]$total_score)))
