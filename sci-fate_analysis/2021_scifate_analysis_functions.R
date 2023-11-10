normal_expr <- function (cds) 
{
  if (!("Size_Factor" %in% names(colData(cds)))) {
    cat("\nSize factor not in cds, calculate size factor...")
    cds = estimate_size_factors(cds)
  }
  b = Matrix::Diagonal(x = 1/cds$Size_Factor)
  return(exprs(cds) %*% b)
}

plot_labels <- function (labels, x, y, cell_size = 0.5, cell_alpha = .5, cols = col_vector, with_labels = T) 
{
  df_cell = data.frame(Cluster = labels, tsne_1 = x, tsne_2 = y)
  cluster_label = df_cell %>% dplyr::group_by(Cluster) %>% 
    dplyr::summarise(mean_1 = median(tsne_1), mean_2 = median(tsne_2))
  if (with_labels) {
    g1 = (ggplot() + geom_point(data = df_cell, aes(x = tsne_1, 
                                                    y = tsne_2, color = Cluster), size = cell_size, stroke = 0, 
                                shape = 16, alpha = cell_alpha) + ggrepel::geom_text_repel(aes(x = cluster_label$mean_1, 
                                                                                      y = cluster_label$mean_2, label = (cluster_label$Cluster))) + 
            xlab("") + ylab("") + 
            # xlim(0,2) + ylim(0,2)+
            theme_void() +
            # coord_fixed() +
            guides(color = guide_legend(title = "", 
                                        override.aes = list(size = 4, alpha = cell_alpha))) + theme(legend.title = element_text(size = 6), 
                                                                                           legend.text = element_text(size = 6), legend.margin = margin(0, 
                                                                                                                                                        -10, 0, 10), legend.key.width = unit(0.15, "in"), 
                                                                                           legend.key.height = unit(0.15, "in"), legend.position = "left"))
  }
  else {
    g1 = (ggplot() + geom_point(data = df_cell, aes(x = tsne_1, 
                                                    y = tsne_2, color = Cluster), size = cell_size, stroke = 0, 
                                shape = 16, alpha = cell_alpha) + xlab("") + ylab("") + 
            # xlim(0,2) + ylim(0,2)+
            theme_void() +
            # coord_fixed() +
            guides(color = guide_legend(title = "", override.aes = list(size = 4, 
                                                                        alpha = 1))) + theme(legend.title = element_text(size = 6), 
                                                                                             legend.text = element_text(size = 6), legend.margin = margin(0, 
                                                                                                                                                          -10, 0, 10), legend.key.width = unit(0.15, "in"), 
                                                                                             legend.key.height = unit(0.15, "in"), legend.position = "left"))
  }
  return(g1)
}

plot_expr <- function (cds, gene, coord_x, coord_y, cell_size = 1, top_v = 3, low_v = -3) 
{
  label_name = gene
  if (gene %in% as.character(rowData(cds)$gene_short_name)) {
    gene_id = as.character((as.data.frame(rowData(cds)) %>% filter(gene_short_name == 
                                                      gene))$gene_id)
    gene_agg_expr = log10(exprs(cds)[gene_id, ]/colData(cds)$Size_Factor + 
                            0.1)
    gene_agg_expr = scale(gene_agg_expr)
    gene_agg_expr[gene_agg_expr > top_v] = top_v
    gene_agg_expr[gene_agg_expr < low_v] = low_v
    colData(cds)$tmp = gene_agg_expr
    
    g1 = (ggplot(data = as.data.frame(colData(cds)), aes(x = coord_x, y = coord_y, 
                                          color = (tmp), alpha = tmp + 1)) 
          + ggrastr::geom_point_rast(size = cell_size, stroke = 0, shape = 16, na.rm = T) 
          + scale_alpha_continuous(range = c(0.5, 1), guide = "none") 
          + theme_void()
          + viridis::scale_color_viridis(option = "viridis", 
                                         name = label_name, na.value = "grey80", 
                                         end = 0.8)
          + guides(alpha = F)
          + theme(legend.title = element_text(size = 6), 
                  legend.text = element_text(size = 6), legend.margin = margin(0, 
                                                                               -10, 0, 10), legend.key.width = unit(0.15, "in"), 
                  legend.key.height = unit(0.15, "in"), legend.position = "left")
          #+ labs(x = "Coordinate 1", y = "Coordinate 2")
    )
    colData(cds)$tmp = NULL
    return(g1)
  }
  else {
    cat("\n gene name: ", gene, " not found")
    return(1)
  }
}


plot_values <- function(values, coord_x, coord_y,  title = "", cell_size = 1) {
  g1 = (ggplot() 
        + ggrastr::geom_point_rast(aes(x = coord_x, y = coord_y, 
                                       color = values), size = cell_size, stroke = 0, shape = 16, na.rm = T) 
        + theme_void()
        + viridis::scale_color_viridis(option = "viridis", 
                                       name = "Values", na.value = "grey80", 
                                       end = 0.8)
        + guides(alpha = F)
        + theme(legend.title = element_text(size = 6), 
                legend.text = element_text(size = 6), legend.margin = margin(0, 
                                                                             -10, 0, 10), legend.key.width = unit(0.15, "in"), 
                legend.key.height = unit(0.15, "in"), legend.position = "left")
        #+ labs(x = "Coordinate 1", y = "Coordinate 2")
        + ggtitle(title)
  )
  return(g1)
  
}

plot_gene_panels <- function (cds, gene_list, coord_x, coord_y, cell_size = 1, top_v = 3, low_v = -3) 
{
  gs = lapply(gene_list, function(x) {
    return(plot_expr(cds, x, coord_x, coord_y, cell_size, top_v, low_v))
  })
  figures = sapply(gs, function(x) {
    typeof(x) != "double"
  })
  gs = gs[figures]
  if(length(gs) > 16) {
    g1 = plot_grid(plotlist = gs, labels = "", ncol = 6, 
                   nrow = 6)
  } else {
    g1 = plot_grid(plotlist = gs, labels = "", ncol = 4, 
                   nrow = 4)
  }
  
  return(g1)
}

plot_gene_module <- function(cds, gene_list, coord_x, coord_y, label_name = "log10(values + 0.1)",cell_size = 1, top_v = 3, low_v = -3) {
  cds_sampled = cds[as.character(gene_list), ]
  gene_sampled = log10(normal_expr(cds_sampled) + 0.1)
  gene_agg_expr = Matrix::colSums(gene_sampled)
  gene_agg_expr = scale(gene_agg_expr)
  gene_agg_expr[gene_agg_expr > top_v] = top_v
  gene_agg_expr[gene_agg_expr < low_v] = low_v
  cds$tmp = gene_agg_expr
  
  g1 = (ggplot(data = colData(cds), aes(x = coord_x, y = coord_y, 
                                        color = (tmp), alpha = tmp + 1)) 
        + ggrastr::geom_point_rast(size = cell_size, stroke = 0, shape = 16, na.rm = T) 
        + viridis::scale_color_viridis(option = "viridis", 
                                       name = label_name, na.value = "grey80", 
                                       end = 0.8) + 
          scale_alpha_continuous(range = c(0.5, 1), guide = "none") 
        + theme_void()
        
        + guides(alpha = F)
        + theme(legend.title = element_text(size = 6), 
                legend.text = element_text(size = 6), legend.margin = margin(0, 
                                                                             -10, 0, 10), legend.key.width = unit(0.15, "in"), 
                legend.key.height = unit(0.15, "in"), legend.position = "none")
        #+ labs(x = "Coordinate 1", y = "Coordinate 2")
  )
  return(g1)
  
}

cal_gene_module <- function(cds, gene_list, top_v = 3, low_v = -3) {
  
  cds_sampled = cds[as.character(gene_list), ]
  gene_sampled = log10(normal_expr(cds_sampled) + 0.1)
  gene_agg_expr = Matrix::colSums(gene_sampled)
  gene_agg_expr = scale(gene_agg_expr)
  gene_agg_expr[gene_agg_expr > top_v] = top_v
  gene_agg_expr[gene_agg_expr < low_v] = low_v
  return(gene_agg_expr)
}

plot_gene_heatmap <- function(cds, top_de_gene_id, sampled_cell_num, ranking = F, cell_rank = "unknown", output_folder = "./", top_v = 3, low_v = -3) {
  if(!file.exists(output_folder)) {
    dir.create(output_folder, recursive = T)
  }
  cds_sampled = cds[as.character(top_de_gene_id), ]
  gene_sampled = log(normal_expr(cds_sampled) + 0.1)
  gene_sampled = t(scale(t(gene_sampled)))
  df_gene = rowData(cds)
  df_gene = df_gene[as.character(top_de_gene_id), ]
  rownames(gene_sampled) = df_gene$gene_short_name
  
  df_cell = colData(cds) %>% select(sample)
  df_cell$cell_rank = cell_rank
  df_cell = df_cell %>% sample_n(sampled_cell_num, replace = F)
  cluster_row = T
  cluster_col = T
  if(ranking == T) {
    cluster_row = F
    cluster_col = F
    df_cell = df_cell %>% arrange(cell_rank)
  }
  gene_sampled = gene_sampled[, as.character(df_cell$sample)]
  
  gene_sampled[gene_sampled > top_v] = top_v
  gene_sampled[gene_sampled < low_v] = low_v
  if(ranking == T) {
    df_anno = df_cell %>% select(cell_rank)
    colnames(df_anno) = "Group"
    rownames(df_anno) = df_cell$sample
    out = pheatmap::pheatmap(gene_sampled, clustering_method="ward.D2", filename = file.path(output_folder, "scheatmap.png"), annotation_col = df_anno,
                             color  = viridis::inferno(100),
                             border_color  = NA,
                             cluster_rows = cluster_row, cluster_cols = cluster_col, 
                             show_rownames = T, show_colnames = F)
  } else {
    out = pheatmap::pheatmap(gene_sampled, clustering_method="ward.D2", filename = file.path(output_folder, "scheatmap.png"), 
                             color  = viridis::inferno(100),
                             border_color  = NA,
                             cluster_rows = cluster_row, cluster_cols = cluster_col, 
                             show_rownames = T, show_colnames = F)
  }
  
  return(out)
  
}

cds_gene_module_identification <- function(cds, top_de_gene_id, sampled_cell_num, plot_coord_x, plot_coord_y, output_folder = "./", cut_tree_rows_n = 20, top_v = 3, low_v = -3) {
  if(!file.exists(output_folder)) {
    dir.create(output_folder, recursive = T)
  }
  
  cds_sampled = cds[as.character(top_de_gene_id), ]
  gene_sampled = log(normal_expr(cds_sampled) + 0.1)
  gene_sampled = t(scale(t(gene_sampled)))
  df_gene = rowData(cds)
  df_gene = df_gene[as.character(top_de_gene_id), ]
  rownames(gene_sampled) = df_gene$gene_short_name
  df_cell = colData(cds) %>% select(sample)
  df_cell = df_cell %>% sample_n(sampled_cell_num, replace = F)
  gene_sampled = gene_sampled[, as.character(df_cell$sample)]
  
  gene_sampled[gene_sampled > top_v] = top_v
  gene_sampled[gene_sampled < low_v] = low_v
  out = pheatmap::pheatmap(gene_sampled, clustering_method="ward.D2", filename = file.path(output_folder, "scheatmap.png"), cutree_rows = cut_tree_rows_n,
                           
                           #clustering_distance_rows = "correlation",
                           #treeheight_row = 10,
                           color  = viridis::inferno(100),
                           border_color  = NA,
                           cluster_rows = T, cluster_cols = T, 
                           show_rownames = F, show_colnames = F)
  
  gene_modules = cutree(out$tree_row, k = cut_tree_rows_n)
  gene_modules = as.data.frame(gene_modules)
  colnames(gene_modules) = "module_id"
  gene_modules$gene_name = rownames(gene_modules)
  gene_modules$gene_id = df_gene$gene_id
  gene_modules = data.frame(gene_modules)
  
  saveRDS(gene_modules, file = file.path(output_folder, "gene_module_list.RDS"))
  
  unique_module = unique(gene_modules$module_id)
  
  plot_list = list()
  id_list = unique_module
  plot_list = lapply(id_list, function(tmp_id) {
    tmp_gene_list = as.character((gene_modules %>% dplyr::filter(module_id == tmp_id))$gene_id)
    cat("\nNumber of genes in module ", tmp_id, " : ", length(tmp_gene_list))
    g1 = plot_gene_module(cds, tmp_gene_list, coord_x = plot_coord_x, coord_y = plot_coord_y, cell_size = 1)
    return(g1)
  })
  
  g1 = plot_grid(plotlist = plot_list, ncol = 6, nrow = 6, labels =  id_list, label_size = 6, label_x = 0, label_y = 1)
  save_plot(g1, filename = file.path(output_folder, "gene_module_exprs.pdf"))
  
  return(list(out, gene_sampled, gene_modules))
}
# plot_expr <- function (cds, gene, coord_x, coord_y, cell_size = 1, top_v = 3, low_v = -3) 
# {
#   label_name = gene
#   if (gene %in% as.character(rowData(cds)$gene_short_name)) {
#     gene_id = as.character((rowData(cds) %>% filter(gene_short_name == 
#                                                       gene))$gene_id)
#     gene_agg_expr = log10(exprs(cds)[gene_id, ]/colData(cds)$Size_Factor + 
#                             0.1)
#     gene_agg_expr = scale(gene_agg_expr)
#     gene_agg_expr[gene_agg_expr > top_v] = top_v
#     gene_agg_expr[gene_agg_expr < low_v] = low_v
#     colData(cds)$tmp = gene_agg_expr
#     
#     g1 = (ggplot(data = colData(cds), aes(x = coord_x, y = coord_y, 
#                                           color = (tmp), alpha = tmp + 1)) 
#           + ggrastr::geom_point_rast(size = cell_size, stroke = 0, shape = 16, na.rm = T) 
#           + scale_alpha_continuous(range = c(0.5, 1), guide = "none") 
#           + theme_void()
#           + viridis::scale_color_viridis(option = "viridis", 
#                                          name = label_name, na.value = "grey80", 
#                                          end = 0.8)
#           + guides(alpha = F)
#           + theme(legend.title = element_text(size = 6), 
#                   legend.text = element_text(size = 6), legend.margin = margin(0, 
#                                                                                -10, 0, 10), legend.key.width = unit(0.15, "in"), 
#                   legend.key.height = unit(0.15, "in"), legend.position = "left")
#           #+ labs(x = "Coordinate 1", y = "Coordinate 2")
#     )
#     colData(cds)$tmp = NULL
#     return(g1)
#   }
#   else {
#     cat("\n gene name: ", gene, " not found")
#     return(1)
#   }
# }

cal_gene_module <- function (cds, gene_list, top_v = 3, low_v = -3) 
{
  cds_sampled = cds[as.character(gene_list), ]
  gene_sampled = log10(normal_expr(cds_sampled) + 0.1)
  gene_agg_expr = Matrix::colSums(gene_sampled)
  gene_agg_expr = scale(gene_agg_expr)
  gene_agg_expr[gene_agg_expr > top_v] = top_v
  gene_agg_expr[gene_agg_expr < low_v] = low_v
  
  result = data.frame("Expression" = gene_agg_expr)
  return(result)
}

cal_gene_module_no_scale <- function (cds, gene_list, top_v = 3, low_v = -3) 
{
  cds_sampled = cds[as.character(gene_list), ]
  gene_sampled = log10(normal_expr(cds_sampled) + 0.1)
  gene_agg_expr = Matrix::colSums(gene_sampled)
  # gene_agg_expr = scale(gene_agg_expr)
  # gene_agg_expr[gene_agg_expr > top_v] = top_v
  # gene_agg_expr[gene_agg_expr < low_v] = low_v
  
  result = data.frame("Expression" = gene_agg_expr)
  return(result)
}
