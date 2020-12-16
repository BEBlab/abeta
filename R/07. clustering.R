#' abetadms__datatable_to_matrix
#'
#' Convert data.table to matrix for heatmap plotting (single AA mutants only).
#'
#' @param input_dt input data.table (required)
#' @param variable_name name of variable to use for heatmap cells (defaut:"fitness")
#'
#' @return a matrix for heamtmap plotting
#' @export
abetadms__datatable_to_matrix <- function(
  input_dt,
  variable_name="nscore_c" 
){
  aa_obj <- Biostrings::AAString("GAVLMIFYWKRHDESTCNQP")
  aa_list <- Biostrings::AMINO_ACID_CODE[strsplit(as.character(aa_obj), NULL)[[1]]]
  aa_list["*"] <- "X"
  
  #Only single AA mutants
  input_df <- as.data.frame(input_dt)
  dms_1aa <- input_df[input_df$Nmut_aa==1,]
  #Absolute position
  # dms_1aa$Pos_abs <- as.numeric(substr(dms_1aa$mut_code, 2, 4))
  #WT sequence
  wt_seq <- unique(dms_1aa[order(dms_1aa[,"Pos_abs"]),c("WT_AA", "Pos_abs")])[,"WT_AA"]
  
  #Construct heatmap matrix
  heat_mat <- matrix(nrow = length(aa_list), ncol = max(dms_1aa$Pos_abs)-min(dms_1aa$Pos_abs)+1)
  rownames(heat_mat) <- names(aa_list)
  colnames(heat_mat) <- min(dms_1aa$Pos_abs):max(dms_1aa$Pos_abs)
  for(aa_pos in min(dms_1aa$Pos_abs):max(dms_1aa$Pos_abs)){
    for(aa_id in names(aa_list)){
      temp_index <- which(dms_1aa$Pos_abs==aa_pos & dms_1aa$Mut==aa_id)
      if(length(temp_index)==1){
        heat_mat[aa_id,as.character(aa_pos)] <- dms_1aa[temp_index,variable_name]
      }
    }
  }
  return(heat_mat)
}



#' abetadms__kmedoids_cluster_singles
#'
#' K-mediods clustering of residue positions based on fitness effects of singles.
#'
#' @param singles_dt data.table with single mutant fitness values (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table


abetadms__kmedoids_cluster_singles <- function(
  singles_dt,
  outpath
){
  
  #Absolute position (singles)
  singles_dt[, Pos_abs := Pos]
  #Mutation code (singles)
  singles_dt[, mut_code := paste0(WT_AA, Pos_abs, Mut)]
  #Mutation code (singles)
  singles_dt[, Nmut_aa := Nham_aa]
  
  #Determine optimimum number of clusters (sample each fitness value from error distribution) - average silhouette width (for K = 1-10)
  set.seed(1)
  pamk_list <- list()
  for(i in 1:100){
    dms_dt_rand <- copy(singles_dt[Nmut_aa==1])
    dms_dt_rand[,nscore_c := rnorm(1, mean = nscore_c, sd = sigma),mut_code]
    heat_mat <- abetadms__datatable_to_matrix(dms_dt_rand)
    d <- dist(t(heat_mat), method = "euclidean") # distance matrix
    pamk_list[[i]] <- fpc::pamk(d,krange=1:10)
  }
  mean_sil_dt <- as.data.table(do.call("rbind", sapply(pamk_list, '[', "crit")))
  names(mean_sil_dt) <- paste0("K=", 1:10)
  #Plot
  plot_dt <- reshape2::melt(mean_sil_dt, measure.vars = 1:10)
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(variable, value)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() + 
    ggplot2::ylab("Average Silhouette width") +
    ggplot2::xlab("Number of clusters")
  ggplot2::ggsave(file=file.path(outpath, 'kmedoids_average_silhouette_width_boxplot.pdf'), width=5, height=5, useDingbats=FALSE)
  
  #Silhouette plot (K=2)
  heat_mat <- abetadms__datatable_to_matrix(singles_dt)
  d <- dist(t(heat_mat), method = "euclidean") # distance matrix
  set.seed(1)
  clust_obj <- fpc::pamk(d,krange=1:10)
  #Plot
  plot_dt <- as.data.table(clust_obj[["pamobject"]][["silinfo"]][["widths"]])
  plot_dt[, Pos := rownames(clust_obj[["pamobject"]][["silinfo"]][["widths"]])]
  plot_dt[, Pos := factor(Pos, levels = Pos)]
  plot_dt[, cluster := factor(cluster)]
  pos_order <- plot_dt[,Pos]
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(Pos, sil_width, fill = cluster)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() +
    ggplot2::ylab("Silhouette width") +
    ggplot2::xlab("Position") +
    ggplot2::scale_fill_grey() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggplot2::ggsave(file=file.path(outpath, 'kmedoids_silhouette_width_barplot.pdf'), width=7, height=5, useDingbats=FALSE)
  
  #Silhouette plot (K=2)
  sil_list <- sapply(sapply(sapply(pamk_list, '[', "pamobject"), '[', "silinfo"), '[', "widths")
  sil_df <- do.call("rbind", sil_list)
  sil_dt <- as.data.table(sil_df)
  sil_dt[, Pos := factor(rownames(sil_df), levels = pos_order)]
  sil_dt[, cluster := factor(cluster)]
  d <- ggplot2::ggplot(sil_dt, ggplot2::aes(Pos, sil_width, fill = cluster)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() + 
    ggplot2::ylab("Silhouette width") +
    ggplot2::xlab("Position") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::facet_grid(cluster~., scales = "free_x")
  ggplot2::ggsave(file=file.path(outpath, 'kmedoids_silhouette_width_boxplots.pdf'), width=7, height=7, useDingbats=FALSE)
  
}





##### run functions
require(data.table)
require(fpc)

dir.create("Clustering")
outpath="Clustering"

#required data:
load("nscore_df.RData")

singles_dt<-setDT(singles)


abetadms__kmedoids_cluster_singles(singles_dt=singles_dt, outpath = outpath)

