library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggupset)

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))


all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = T, sep='\t')
all_genes = all_genes[!(all_genes$tissue=="MCA" & all_genes$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]
all_genes = all_genes[!(all_genes$tissue=="MCA" & all_genes$project %in% c("AngII", "Salt-sensitive")),]
all_genes$comp = paste(all_genes$strain, all_genes$treatment,  sep="-")

deg_count = table(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), ]$gene_name)
deg_tissue_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), c("gene_name", "tissue")])$gene_name)
deg_strain_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), c("gene_name", "strain")])$gene_name)

analyze_deg_overlap_table <- function(data, group_var = "tissue", 
                                      strain_filter = c("C57BL/6", "SS", "SHR"), 
                                      pval_thresh = 0.05, logfc_thresh = 0.5, 
                                      n_sim = 1000, seed = 42) {
  set.seed(seed)
  group_var <- rlang::sym(group_var)
  
  # Filter real DEGs
  deg_filtered <- data %>%
    filter(p_val_adj < pval_thresh,
           abs(avg_log2FC) > logfc_thresh,
           strain %in% strain_filter)
  
  # Observed: how many groups each gene is DEG in
  deg_gene_group <- deg_filtered %>%
    distinct(gene_name, !!group_var)
  
  gene_group_table_obs <- table(deg_gene_group$gene_name)
  freq_obs <- as.vector(table(gene_group_table_obs))  # how many genes appear in 1, 2, 3... groups
  names(freq_obs) <- names(table(gene_group_table_obs))
  
  # Prepare for simulation
  # all_genes <- unique(data$gene_name)
  all_genes <- unique(data$gene_name)
  
  # DEG counts per group
  deg_count_per_group <- deg_gene_group %>%
    dplyr::count(!!group_var, name = "n_deg")
  
  # Simulate
  sim_freq_list <- replicate(n_sim, {
    sim_deg_gene_group <- deg_count_per_group %>%
      mutate(gene_name = map(n_deg, ~ sample(all_genes, .x))) %>%
      unnest(gene_name)
    
    gene_group_table_sim <- table(sim_deg_gene_group$gene_name)
    freq_sim <- table(gene_group_table_sim)
    
    # Ensure all group sizes are present
    out <- rep(0, max(length(freq_obs), as.numeric(names(freq_sim)) %>% max()))
    names(out) <- as.character(1:length(out))
    out[names(freq_sim)] <- as.numeric(freq_sim)
    out
  }, simplify = "matrix")
  
  sim_freq_mean <- rowMeans(sim_freq_list)
  sim_freq_sd <- apply(sim_freq_list, 1, sd)
  z_scores <- (freq_obs - sim_freq_mean) / sim_freq_sd
  
  result <- list(
    observed_counts = freq_obs,
    sim_mean = sim_freq_mean,
    sim_sd = sim_freq_sd,
    z_score = z_scores,
    sim_matrix = sim_freq_list
  )
  
  return(result)
}



result_tissue <- analyze_deg_overlap_table(all_genes, group_var = "tissue")
result_strain <- analyze_deg_overlap_table(all_genes, group_var = "strain")

result_tissue$z_score
result_strain$z_score


result = result_tissue
observed_kplus <- sum(result$observed_counts[names(result$observed_counts) >= 2])
sim_kplus <- colSums(result$sim_matrix[as.numeric(rownames(result$sim_matrix)) >= 2, ])
result$p_empirical <- mean(sim_kplus >= observed_kplus)
result$p_text <- paste0("p-value = ", signif(result$p_empirical, 2), 
                 "\n(1,000 perms)")
result_tissue=result


result = result_strain
observed_kplus <- sum(result$observed_counts[names(result$observed_counts) >= 2])
sim_kplus <- colSums(result$sim_matrix[as.numeric(rownames(result$sim_matrix)) >= 2, ])
result$p_empirical <- mean(sim_kplus >= observed_kplus)
result$p_text <- paste0("p-value = ", signif(result$p_empirical, 2), 
                        "\n(1,000 perms)")
result_strain=result



result = result_tissue; group_type="tissues"
plot_df <- tibble(
  GroupCount = as.numeric(names(result$observed_counts)),
  Observed = as.numeric(result$observed_counts),
  Expected = result$sim_mean[1:length(result$observed_counts)],
  SD = result$sim_sd[1:length(result$observed_counts)]
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
# write.table(plot_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.count.tissue.obs_exp.tsv", sep = "\t", col.names = T, row.names = F)
plot_long <- plot_df %>%
  dplyr::select(GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")
x_pos <- max(plot_df$GroupCount) + 0.5
y_pos <- max(plot_df$Expected) * 0.95
p_val <- result$p_text
p1<-ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  # geom_line(aes(y = Expected), color = "gray40", size = 1) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray40") +
  geom_line(data = plot_long, aes(y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray40")) +
  annotate("text", x = x_pos, y = y_pos, label = p_val,
           hjust = 1, vjust = 1, size = 4.2) +
  labs(
    x = paste("Number of shared", group_type),
    y = "Number of genes"
  )

result = result_strain; group_type="strains"
plot_df <- tibble(
  GroupCount = as.numeric(names(result$observed_counts)),
  Observed = as.numeric(result$observed_counts),
  Expected = result$sim_mean[1:length(result$observed_counts)],
  SD = result$sim_sd[1:length(result$observed_counts)]
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
# write.table(plot_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.count.strain.obs_exp.tsv", sep = "\t", col.names = T, row.names = F)
plot_long <- plot_df %>%
  dplyr::select(GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")
x_pos <- max(plot_df$GroupCount) + 0.5
y_pos <- max(plot_df$Expected) * 0.95
p_val <- result$p_text
p2<-ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  # geom_line(aes(y = Expected), color = "gray40", size = 1) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray40") +
  geom_line(data = plot_long, aes(y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray40")) +
  annotate("text", x = x_pos, y = y_pos, label = p_val,
           hjust = 1, vjust = 1, size = 4.2) +
  labs(
    x = paste("Number of shared", group_type),
    y = "Number of genes"
  )

p1+p2

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.obs_exp.png", width=793/96, height=220/96, dpi=300)


result_tissue_hypt <- analyze_deg_overlap_table(all_genes, group_var = "tissue")
result_strain_hypt <- analyze_deg_overlap_table(all_genes, group_var = "strain")
result_tissue_hypt$z_score
result_strain_hypt$z_score

result_tissue_norm <- analyze_deg_overlap_table(all_genes, group_var = "tissue", strain_filter = c("SD", "WKY"))
result_strain_norm <- analyze_deg_overlap_table(all_genes, group_var = "strain", strain_filter = c("SD", "WKY"))
result_tissue_norm$z_score
result_strain_norm$z_score

z_hypt_tissue <- result_tissue_hypt$z_score
z_norm_tissue <- result_tissue_norm$z_score
z_hypt_strain <- result_strain_hypt$z_score
z_norm_strain <- result_strain_norm$z_score

df_zscore <- bind_rows(
  tibble(Group = as.integer(names(z_hypt_tissue)),
         Z = as.numeric(z_hypt_tissue),
         Context = "Hypertensive",
         Type = "Tissue"),
  tibble(Group = as.integer(names(z_norm_tissue)),
         Z = as.numeric(z_norm_tissue),
         Context = "Normotensive",
         Type = "Tissue"),
  tibble(Group = as.integer(names(z_hypt_strain)),
         Z = as.numeric(z_hypt_strain),
         Context = "Hypertensive",
         Type = "Strain"),
  tibble(Group = as.integer(names(z_norm_strain)),
         Z = as.numeric(z_norm_strain),
         Context = "Normotensive",
         Type = "Strain")
) %>%
  complete(Group, Type, Context, fill = list(Z = NA))

df_zscore$Z[is.infinite(df_zscore$Z)] <- max(df_zscore$Z[is.finite(df_zscore$Z)], na.rm = TRUE) + 10
df_zscore <- df_zscore[!(df_zscore$Type=="Strain" & df_zscore$Group >3), ]

ggplot(df_zscore, aes(x = factor(Group), y = Z, fill = Context)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Type, scales = "free_x", nrow = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(values = c("Hypertensive" = "#D73027", "Normotensive" = "#4575B4")) +
  labs(
    x = "Number of shared tissues or strains",
    y = "Z-score of DEG sharing\n(obs/exp)",
    fill = "Strain"
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.obs_exp.hypt_norm.png", width=480/96, height=297/96, dpi=300)






analyze_deg_overlap_within_tissue <- function(data,
                                              strain_filter = c("C57BL/6", "SS", "SHR"), 
                                              pval_thresh = 0.05,
                                              logfc_thresh = 0.5,
                                              n_sim = 1000,
                                              seed = 42) {
  set.seed(seed)
  
  all_genes <- unique(data$gene_name)
  tissue_list <- unique(data$tissue)
  
  # Prepare storage
  result_list <- list()
  sim_matrix_list <- list()
  
  for (tissue in tissue_list) {
    df <- data %>%
      filter(tissue == !!tissue,
             p_val_adj < pval_thresh,
             abs(avg_log2FC) > logfc_thresh,
             strain %in% strain_filter) # %>%
      distinct(gene_name, cell_type)
    
    gene_cell_table <- table(df$gene_name)
    obs_freq <- table(gene_cell_table)
    
    deg_counts <- df %>%
      group_by(cell_type) %>%
      summarise(n_deg = n(), .groups = "drop")
    
    sim_matrix <- replicate(n_sim, {
      sim_df <- deg_counts %>%
        mutate(gene_name = map(n_deg, ~ sample(all_genes, .x))) %>%
        unnest(gene_name)
      
      gene_cell_table_sim <- table(sim_df$gene_name)
      sim_freq <- table(gene_cell_table_sim)
      
      out <- rep(0, max(c(as.integer(names(obs_freq)), as.integer(names(sim_freq)))))
      names(out) <- as.character(seq_along(out))
      out[names(sim_freq)] <- sim_freq
      out
    }, simplify = "matrix")
    
    sim_mean <- rowMeans(sim_matrix)
    sim_sd <- apply(sim_matrix, 1, sd)
    obs_counts <- as.numeric(obs_freq)
    names(obs_counts) <- names(obs_freq)
    
    n_bins <- max(length(obs_counts), length(sim_mean))
    obs_counts <- `length<-`(obs_counts, n_bins)
    sim_mean <- `length<-`(sim_mean, n_bins)
    sim_sd <- `length<-`(sim_sd, n_bins)
    
    z_scores <- (obs_counts - sim_mean) / sim_sd
    
    result_df <- tibble(
      tissue = tissue,
      n_celltypes_shared = seq_along(z_scores),
      z_score = z_scores,
      observed = obs_counts,
      expected = sim_mean,
      sd = sim_sd
    )
    
    result_list[[tissue]] <- result_df
    sim_matrix_list[[tissue]] <- sim_matrix
  }
  
  return(list(
    result = bind_rows(result_list),
    sim_matrix = sim_matrix_list
  ))
}

out <- analyze_deg_overlap_within_tissue(all_genes)


sim_matrix_list <- out$sim_matrix
result_within_tissue <- out$result

pval_table

pval_positions <- plot_df %>%
  group_by(Tissue) %>%
  summarise(
    x_pos = max(GroupCount, na.rm = TRUE)) %>%
  left_join(pval_table, by = c("Tissue" = "tissue")) %>%
  mutate(
    label = paste0("p = ", formatC(empirical_p, format = "e", digits = 2)),
    y_pos = max(plot_df$Expected) * 0.9  # position the label above x=2
  )

result <- out$result
plot_df <- tibble(
  Tissue = result$tissue,
  GroupCount = as.numeric(result$n_celltypes_shared),
  Observed = as.numeric(result$observed),
  Expected = as.numeric(result$expected),
  SD = as.numeric(result$sd)
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
plot_long <- plot_df %>%
  dplyr::select(Tissue, GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")

ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  # geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray80") +
  # geom_line(data = plot_long, aes(x = as.numeric(GroupCount), y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray80")) +
  geom_text(aes(y = Observed, label = Observed),
            vjust = -0.5, size = 3, color = "black") +
  # geom_text(data = pval_positions,
            # aes(x = x_pos, y = y_pos, label = label),
            # inherit.aes = FALSE,
            # size = 3.5, hjust = 1) +
  labs(
    x = "Number of shared cell types",
    y = "Number of genes"
  ) +
  facet_grid(~Tissue, scales = "free_x", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.tissue_wise.w_strain.obs_exp.png", width=1200/96, height=220/96, dpi=300)



result$label = ifelse(is.infinite(result$z_score),
                      ifelse(result$z_score > 0, "Inf", "-Inf"),
                      round(result$z_score, 1))
result$y = ifelse(is.infinite(result$z_score),
                  ifelse(result$z_score > 0, max(result$z_score[is.finite(result$z_score)], na.rm = TRUE) + 100,
                         min(result$z_score[is.finite(result$z_score)], na.rm = TRUE) - 10),
                  result$z_score)
ggplot(result, aes(x = factor(n_celltypes_shared), y = z_score)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "lightgrey", fill = "lightgrey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylim(-90, 1400) +
  labs(
    x = "Number of shared cell types",
    y = "Z-score of DEG sharing\n(obs/exp)",
    fill = "Strain"
  ) +
  geom_text(data = result,
            aes(label = label, y = y),
            vjust = ifelse(result$z_score < 0, 1.1, -0.5),
            size = 3) +
  facet_grid(~tissue, scales = "free_x", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.tissue_wise.w_strain.zscore.barplot.png", width=1100/96, height=260/96, dpi=300)




# tissue observed_shared_deg_ge2 empirical_p
# <chr>                    <dbl>       <dbl>
#   1 HYP                        988           0
# 2 LK                         314           0
# 3 LV                         229           0
# 4 MCA                         10           0
# 5 MSA                         15           0
# 

















################################################################################
### abundance matched 

deg_pseudo = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/pseudo.DEG.all.out", sep='\t', header=T)
colnames(deg_pseudo)[12] = "strain"
deg_pseudo$cell_type = factor(deg_pseudo$cell_type, levels = cell_order)
deg_pseudo$tissue = factor(deg_pseudo$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_pseudo$treatment = factor(deg_pseudo$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_pseudo$strain = factor(deg_pseudo$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

expr_all = unique(rbind(data.frame(unname(deg_pseudo[, c("pct.1", "avg_expr.1", "gene_name", "cell_type", "project", "strain", "tissue", "control")])),
                        data.frame(unname(deg_pseudo[, c("pct.2", "avg_expr.2", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")]))))
names(expr_all) = c("pct", "avg_expr", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")
expr_all$cell_type = factor(expr_all$cell_type, levels = cell_order)
expr_all$tissue = factor(expr_all$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
expr_all$treatment = factor(expr_all$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
expr_all$strain = factor(expr_all$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))



gene_abundance_tile <- expr_all %>%
  group_by(gene_name, cell_type, tissue) %>%
  summarise(pct_avg = mean(pct),
            expr_avg = mean(avg_expr),
            abundance_score = pct_avg * expr_avg,
            .groups = "drop") %>%
  group_by(cell_type, tissue) %>%
  mutate(abundance_bin = ntile(abundance_score, 3)) %>%
  ungroup()

gene_abundance_threshold <- expr_all %>%
  group_by(gene_name, cell_type, tissue) %>%
  summarise(pct_avg = mean(pct),
            expr_avg = mean(avg_expr),
            abundance_score = pct_avg * expr_avg,
            .groups = "drop") %>%
  group_by(cell_type, tissue) %>%
  mutate(abundance_bin = case_when(
    pct_avg <= 0.3 ~ "Low",
    pct_avg <= 0.6 ~ "Medium",
    pct_avg > 0.6  ~ "High"
  )) %>%
  ungroup()

gene_abundance = gene_abundance_tile

all_genes_by_bin_cell_type <- gene_abundance %>%
  distinct(gene_name, tissue, cell_type, abundance_bin) %>%
  mutate(key = paste(tissue, cell_type, abundance_bin, sep = "___")) %>%
  group_by(key) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  deframe()

all_genes_by_bin_tissue <- gene_abundance %>%
  distinct(gene_name, tissue, abundance_bin) %>%
  mutate(key = paste(tissue, abundance_bin, sep = "___")) %>%
  group_by(key) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  deframe()

all_gene_by_bin = all_genes_by_bin_cell_type

all_genes <- all_genes %>%
  left_join(gene_abundance %>% dplyr::select(gene_name, tissue, cell_type, abundance_bin), 
            by = c("gene_name", "tissue", "cell_type"))

analyze_deg_overlap_matched_abundance <- function(data, group_var = "tissue", 
                                                  strain_filter = c("C57BL/6", "SS", "SHR"), 
                                                  pval_thresh = 0.05, logfc_thresh = 0.5, 
                                                  n_sim = 1000, seed = 42) {
  set.seed(seed)
  group_var_sym <- rlang::sym(group_var)
  
  # Step 1: Filter DEGs
  deg_filtered <- data %>%
    filter(p_val_adj < pval_thresh,
           abs(avg_log2FC) > logfc_thresh,
           strain %in% strain_filter)
  
  # Step 2: Observed count of DEG per gene across groups
  deg_gene_group <- deg_filtered %>%
    distinct(gene_name, !!group_var_sym)
  
  gene_group_table_obs <- table(deg_gene_group$gene_name)
  freq_obs <- table(gene_group_table_obs)
  freq_obs_vec <- rep(0, max(as.integer(names(freq_obs))))
  names(freq_obs_vec) <- as.character(seq_along(freq_obs_vec))
  freq_obs_vec[names(freq_obs)] <- as.numeric(freq_obs)
  
  # Normalize observed to relative proportions
  prop_obs_vec <- freq_obs_vec / sum(freq_obs_vec)
  
  # Step 3: Count DEGs by group × abundance_bin
  deg_count_per_group_bin <- deg_filtered %>%
    distinct(gene_name, !!group_var_sym, tissue, cell_type, abundance_bin) %>%
    dplyr::count(!!group_var_sym, cell_type, tissue, abundance_bin)
  
  # Step 4: Simulate null distribution
  sim_matrix_prop <- replicate(n_sim, {
    sim_deg_gene_group <- deg_count_per_group_bin %>%
      mutate(gene_name = pmap(list(n, tissue, cell_type, abundance_bin), function(n_genes, tiss, ct, bin) {
        key <- paste(tiss, ct, bin, sep = "___")
        gene_pool <- all_gene_by_bin[[key]]
        sample(gene_pool, n_genes)
      })) %>%
      unnest(gene_name)
    
    sim_deg_gene_group <- sim_deg_gene_group %>% distinct(gene_name, !!group_var_sym)
    gene_group_table_sim <- table(sim_deg_gene_group$gene_name)
    freq_sim <- table(gene_group_table_sim)
    
    sim_vec <- rep(0, length(freq_obs_vec))
    names(sim_vec) <- names(freq_obs_vec)
    sim_vec[names(freq_sim)] <- as.numeric(freq_sim)
    sim_vec / sum(sim_vec)  # Normalize to proportion
  }, simplify = "matrix")
  
  # Step 5: Compute mean, sd, z-score from proportions
  sim_mean <- rowMeans(sim_matrix_prop)
  sim_sd <- apply(sim_matrix_prop, 1, sd)
  z_scores <- (prop_obs_vec - sim_mean) / sim_sd
  
  list(
    observed_proportion = prop_obs_vec,
    sim_mean = sim_mean,
    sim_sd = sim_sd,
    z_score = z_scores,
    sim_matrix = sim_matrix_prop
  )
}







result_tissue <- analyze_deg_overlap_matched_abundance(all_genes, group_var = "tissue")
result_strain <- analyze_deg_overlap_matched_abundance(all_genes, group_var = "strain")

result_tissue$z_score
result_strain$z_score


result = result_tissue
observed_kplus <- sum(result$observed_proportion[names(result$observed_proportion) >= 3])
sim_kplus <- colSums(result$sim_matrix[as.numeric(rownames(result$sim_matrix)) >= 3, ])
result$p_empirical <- mean(sim_kplus >= observed_kplus)
result$p_text <- paste0("p-value = ", signif(result$p_empirical, 2), 
                        "\n(≥3 groups,\n1,000 perms)")
result_tissue=result


result = result_strain
observed_kplus <- sum(result$observed_proportion[names(result$observed_proportion) >= 3])
sim_kplus <- result$sim_matrix[as.numeric(rownames(result$sim_matrix)) >= 3, ]
result$p_empirical <- mean(sim_kplus >= observed_kplus)
result$p_text <- paste0("p-value = ", signif(result$p_empirical, 2), 
                        "\n(≥3 groups,\n1,000 perms)")
result_strain=result



result = result_tissue; group_type="tissues"
plot_df <- tibble(
  GroupCount = as.numeric(names(result$observed_proportion)),
  Observed = as.numeric(result$observed_proportion),
  Expected = result$sim_mean[1:length(result$observed_proportion)],
  SD = result$sim_sd[1:length(result$observed_proportion)]
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
write.table(plot_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.count.tissue.obs_exp.expr_match.tsv", sep = "\t", col.names = T, row.names = F)
plot_long <- plot_df %>%
  dplyr::select(GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")
x_pos <- max(plot_df$GroupCount) + 0.5
y_pos <- max(plot_df$Expected) * 0.95
p_val <- result$p_text
p1<-ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  # geom_line(aes(y = Expected), color = "gray40", size = 1) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray40") +
  geom_line(data = plot_long, aes(y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray40")) +
  # annotate("text", x = x_pos, y = y_pos, label = p_val,
  #          hjust = 1, vjust = 1, size = 4.2) +
  labs(
    x = paste("Number of shared", group_type),
    y = "Proportion of genes"
  )

result = result_strain; group_type="strains"
plot_df <- tibble(
  GroupCount = as.numeric(names(result$observed_proportion)),
  Observed = as.numeric(result$observed_proportion),
  Expected = result$sim_mean[1:length(result$observed_proportion)],
  SD = result$sim_sd[1:length(result$observed_proportion)]
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
write.table(plot_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.count.strain.obs_exp.expr_match.tsv", sep = "\t", col.names = T, row.names = F)
plot_long <- plot_df %>%
  dplyr::select(GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")
x_pos <- max(plot_df$GroupCount) + 0.5
y_pos <- max(plot_df$Expected) * 0.95
p_val <- result$p_text
p2<-ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  # geom_line(aes(y = Expected), color = "gray40", size = 1) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray40") +
  geom_line(data = plot_long, aes(y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray40")) +
  annotate("text", x = x_pos, y = y_pos, label = p_val,
           hjust = 1, vjust = 1, size = 4) +
  scale_x_continuous(breaks=c(1,2,3)) +
  labs(
    x = paste("Number of shared", group_type),
    y = "Proportion of genes"
  )

p1+p2

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.obs_exp.expr_match.png", width=793/96, height=220/96, dpi=300)




### group >= 3 only
result = result_tissue; group_type="tissues"
plot_df <- tibble(
  GroupCount = as.numeric(names(result$observed_proportion)),
  Observed = as.numeric(result$observed_proportion),
  Expected = result$sim_mean[1:length(result$observed_proportion)],
  SD = result$sim_sd[1:length(result$observed_proportion)]
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  ) %>%
  filter(GroupCount>2) 
plot_long <- plot_df %>%
  dplyr::select(GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")
x_pos <- max(plot_df$GroupCount) + 0.5
y_pos <- max(plot_df$Observed)
p_val <- result$p_text
ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  # geom_line(aes(y = Expected), color = "gray40", size = 1) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray40") +
  geom_line(data = plot_long, aes(y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray40")) +
  theme(legend.position = "None") +
  scale_x_continuous(breaks=c(3, 4, 5)) +
  annotate("text", x = x_pos, y = y_pos, label = p_val,
           hjust = 1, vjust = 1, size = 4) +
  labs(
    x = "",
    y = ""
  )
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.obs_exp.expr_match.g_3.png", width=245/96, height=180/96, dpi=300)










result_tissue_hypt <- analyze_deg_overlap_matched_abundance(all_genes, group_var = "tissue")
result_strain_hypt <- analyze_deg_overlap_matched_abundance(all_genes, group_var = "strain")
result_tissue_hypt$z_score
result_strain_hypt$z_score

result_tissue_norm <- analyze_deg_overlap_matched_abundance(all_genes, group_var = "tissue", strain_filter = c("SD", "WKY"))
result_strain_norm <- analyze_deg_overlap_matched_abundance(all_genes, group_var = "strain", strain_filter = c("SD", "WKY"))
result_tissue_norm$z_score
result_strain_norm$z_score

z_hypt_tissue <- result_tissue_hypt$z_score
z_norm_tissue <- result_tissue_norm$z_score
z_hypt_strain <- result_strain_hypt$z_score
z_norm_strain <- result_strain_norm$z_score

df_zscore <- bind_rows(
  tibble(Group = as.integer(names(z_hypt_tissue)),
         Z = as.numeric(z_hypt_tissue),
         Context = "Hypertensive",
         Type = "Tissue"),
  tibble(Group = as.integer(names(z_norm_tissue)),
         Z = as.numeric(z_norm_tissue),
         Context = "Normotensive",
         Type = "Tissue"),
  tibble(Group = as.integer(names(z_hypt_strain)),
         Z = as.numeric(z_hypt_strain),
         Context = "Hypertensive",
         Type = "Strain"),
  tibble(Group = as.integer(names(z_norm_strain)),
         Z = as.numeric(z_norm_strain),
         Context = "Normotensive",
         Type = "Strain")
) %>%
  complete(Group, Type, Context, fill = list(Z = NA))

df_zscore$Z[is.infinite(df_zscore$Z)] <- max(df_zscore$Z[is.finite(df_zscore$Z)], na.rm = TRUE) + 10
df_zscore <- df_zscore[!(df_zscore$Type=="Strain" & df_zscore$Group >3), ]

ggplot(df_zscore, aes(x = factor(Group), y = Z, fill = Context)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Type, scales = "free_x", nrow = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(values = c("Hypertensive" = "#D73027", "Normotensive" = "#4575B4")) +
  labs(
    x = "Number of shared tissues or strains",
    y = "Z-score of DEG sharing\n(obs/exp)",
    fill = "Strain"
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.obs_exp.hypt_norm.expr_match.png", width=480/96, height=297/96, dpi=300)







analyze_deg_overlap_within_tissue_abundance <- function(data,
                                                        strain_filter = c("C57BL/6", "SS", "SHR"), 
                                                        pval_thresh = 0.05,
                                                        logfc_thresh = 0.5,
                                                        n_sim = 1000,
                                                        seed = 42) {
  set.seed(seed)
  tissue_list <- unique(data$tissue)
  result_list <- list()
  sim_matrix_list <- list()
  
  for (tissue in tissue_list) {
    deg_filtered <- data %>%
      filter(tissue == !!tissue,
             p_val_adj < pval_thresh,
             abs(avg_log2FC) > logfc_thresh,
             strain %in% strain_filter)
    
    deg_gene_group <- deg_filtered %>%
      distinct(gene_name, cell_type)
    
    gene_group_table_obs <- table(deg_gene_group$gene_name)
    freq_obs <- table(gene_group_table_obs)
    freq_obs_vec <- rep(0, max(as.integer(names(freq_obs))))
    names(freq_obs_vec) <- as.character(seq_along(freq_obs_vec))
    freq_obs_vec[names(freq_obs)] <- as.numeric(freq_obs)
    prop_obs_vec <- freq_obs_vec / sum(freq_obs_vec)
    
    deg_count_per_bin <- deg_filtered %>%
      distinct(gene_name, cell_type, tissue, abundance_bin) %>%
      dplyr::count(cell_type, tissue, abundance_bin)
    
    sim_matrix_prop <- replicate(n_sim, {
      sim_deg_gene_group <- deg_count_per_bin %>%
        mutate(gene_name = pmap(list(n, tissue, cell_type, abundance_bin), function(n_genes, tiss, ct, bin) {
          key <- paste(tiss, ct, bin, sep = "___")
          gene_pool <- all_gene_by_bin[[key]]
          sample(gene_pool, n_genes, replace = TRUE)
        })) %>%
        unnest(gene_name)
      
      sim_deg_gene_group <- sim_deg_gene_group %>% distinct(gene_name, cell_type)
      gene_group_table_sim <- table(sim_deg_gene_group$gene_name)
      freq_sim <- table(gene_group_table_sim)
      
      sim_vec <- rep(0, length(freq_obs_vec))
      names(sim_vec) <- names(freq_obs_vec)
      sim_vec[names(freq_sim)] <- as.numeric(freq_sim)
      sim_vec / sum(sim_vec)
    }, simplify = "matrix")
    
    sim_mean <- rowMeans(sim_matrix_prop)
    sim_sd <- apply(sim_matrix_prop, 1, sd)
    z_scores <- (prop_obs_vec - sim_mean) / sim_sd
    
    result_df <- tibble(
      tissue = tissue,
      n_celltypes_shared = seq_along(z_scores),
      z_score = z_scores,
      observed = prop_obs_vec,
      expected = sim_mean,
      sd = sim_sd
    )
    
    result_list[[tissue]] <- result_df
    sim_matrix_list[[tissue]] <- sim_matrix_prop
  }
  
  return(list(
    result = bind_rows(result_list),
    sim_matrix = sim_matrix_list
  ))
}

out <- analyze_deg_overlap_within_tissue_abundance(all_genes)


sim_matrix_list <- out$sim_matrix
result_within_tissue <- out$result

pval_table <- map_df(names(sim_matrix_list), function(tissue) {
  
  # Get simulation matrix and observed values for this tissue
  sim_mat <- sim_matrix_list[[tissue]]
  obs_vals <- result_within_tissue %>%
    filter(!(is.na(observed))) %>%
    filter(tissue == !!tissue, n_celltypes_shared >= 2) %>%
    pull(observed)
  
  observed_total <- sum(obs_vals)
  
  sim_totals <- colSums(sim_mat[as.integer(rownames(sim_mat)) >= 2, , drop = FALSE])
  
  p_emp <- mean(sim_totals >= observed_total)
  
  tibble(
    tissue = tissue,
    observed_shared_deg_ge2 = observed_total,
    empirical_p = p_emp
  )
})

pval_table

pval_positions <- plot_df %>%
  group_by(Tissue) %>%
  summarise(
    x_pos = max(GroupCount, na.rm = TRUE)) %>%
  left_join(pval_table, by = c("Tissue" = "tissue")) %>%
  mutate(
    label = paste0("p = ", formatC(empirical_p, format = "e", digits = 2)),
    y_pos = max(plot_df$Expected) * 0.9  # position the label above x=2
  )

result <- out$result
plot_df <- tibble(
  Tissue = result$tissue,
  GroupCount = as.numeric(result$n_celltypes_shared),
  Observed = as.numeric(result$observed),
  Expected = as.numeric(result$expected),
  SD = as.numeric(result$sd)
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
plot_long <- plot_df %>%
  dplyr::select(Tissue, GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")

ggplot(plot_df, aes(x = factor(GroupCount))) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray80") +
  geom_line(data = plot_long, aes(x = as.numeric(GroupCount), y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray80")) +
  geom_text(aes(y = Observed, label = Observed),
            vjust = -0.5, size = 3, color = "black") +
  # geom_text(data = pval_positions,
  #           aes(x = x_pos, y = y_pos, label = label),
  #           inherit.aes = FALSE,
  #           size = 3.5, hjust = 1) +
  labs(
    x = "Number of shared cell types",
    y = "Number of genes"
  ) +
  facet_grid(~Tissue, scales = "free_x", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.tissue_wise.w_strain.obs_exp.expr_match.png", width=1200/96, height=220/96, dpi=300)



result$label = ifelse(is.infinite(result$z_score),
                      ifelse(result$z_score > 0, "Inf", "-Inf"),
                      round(result$z_score, 1))
result$y = ifelse(is.infinite(result$z_score),
                  ifelse(result$z_score > 0, max(result$z_score[is.finite(result$z_score)], na.rm = TRUE) + 100,
                         min(result$z_score[is.finite(result$z_score)], na.rm = TRUE) - 10),
                  result$z_score)
ggplot(result, aes(x = factor(n_celltypes_shared), y = z_score)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "lightgrey", fill = "lightgrey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylim(-50, 1400) +
  labs(
    x = "Number of shared cell types",
    y = "Z-score of DEG sharing\n(obs/exp)",
    fill = "Strain"
  ) +
  geom_text(data = result,
            aes(label = label, y = y),
            vjust = ifelse(result$z_score < 0, 1.1, -0.5),
            size = 3) +
  facet_grid(~tissue, scales = "free_x", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.tissue_wise.w_strain.zscore.barplot.expr_match.png", width=1100/96, height=260/96, dpi=300)





















analyze_deg_overlap_within_tissue_abundance <- function(data,
                                                        strain_filter = c("C57BL/6", "SS", "SHR"), 
                                                        pval_thresh = 0.05,
                                                        logfc_thresh = 0.5,
                                                        n_sim = 1000,
                                                        seed = 42) {
  set.seed(seed)
  tissue_list <- unique(data$tissue)
  result_list <- list()
  sim_matrix_list <- list()
  
  for (tissue in tissue_list) {
    deg_filtered <- data %>%
      filter(tissue == !!tissue,
             p_val_adj < pval_thresh,
             abs(avg_log2FC) > logfc_thresh,
             strain %in% strain_filter)
    
    deg_gene_group <- deg_filtered %>%
      distinct(gene_name, cell_type)
    
    gene_group_table_obs <- table(deg_gene_group$gene_name)
    freq_obs <- table(gene_group_table_obs)
    
    # Collapse into 2 groups: 1 vs. 2+
    group_labels <- ifelse(as.integer(names(freq_obs)) == 1, "1", "2+")
    freq_obs_collapsed <- tapply(as.numeric(freq_obs), group_labels, sum)
    prop_obs_vec <- freq_obs_collapsed / sum(freq_obs_collapsed)
    
    deg_count_per_bin <- deg_filtered %>%
      distinct(gene_name, cell_type, tissue, abundance_bin) %>%
      dplyr::count(cell_type, tissue, abundance_bin)
    
    sim_matrix_prop <- replicate(n_sim, {
      sim_deg_gene_group <- deg_count_per_bin %>%
        mutate(gene_name = pmap(list(n, tissue, cell_type, abundance_bin), function(n_genes, tiss, ct, bin) {
          key <- paste(tiss, ct, bin, sep = "___")
          gene_pool <- all_gene_by_bin[[key]]
          sample(gene_pool, n_genes, replace = TRUE)
        })) %>%
        unnest(gene_name)
      
      sim_deg_gene_group <- sim_deg_gene_group %>% distinct(gene_name, cell_type)
      gene_group_table_sim <- table(sim_deg_gene_group$gene_name)
      freq_obs_sim <- table(gene_group_table_sim)
      
      group_labels_sim <- ifelse(as.integer(names(freq_obs_sim)) == 1, "1", "2+")
      freq_sim_collapsed <- tapply(as.numeric(freq_obs_sim), group_labels_sim, sum)
      prop_sim_vec <- freq_sim_collapsed / sum(freq_sim_collapsed)
      
      # Ensure the same order and length
      out_vec <- prop_obs_vec
      out_vec[names(prop_sim_vec)] <- prop_sim_vec
      out_vec
    }, simplify = "matrix")
    
    sim_mean <- rowMeans(sim_matrix_prop)
    sim_sd <- apply(sim_matrix_prop, 1, sd)
    z_scores <- (prop_obs_vec - sim_mean) / sim_sd
    
    result_df <- tibble(
      tissue = tissue,
      group = names(prop_obs_vec),
      z_score = z_scores,
      observed = prop_obs_vec,
      expected = sim_mean,
      sd = sim_sd
    )
    
    result_list[[tissue]] <- result_df
    sim_matrix_list[[tissue]] <- sim_matrix_prop
  }
  
  return(list(
    result = bind_rows(result_list),
    sim_matrix = sim_matrix_list
  ))
}

out <- analyze_deg_overlap_within_tissue_abundance(all_genes)


sim_matrix_list <- out$sim_matrix
result_within_tissue <- out$result

pval_table <- map_df(names(sim_matrix_list), function(tissue) {
  
  # Get simulation matrix and observed values for this tissue
  sim_mat <- sim_matrix_list[[tissue]]
  obs_vals <- result_within_tissue %>%
    filter(!(is.na(observed))) %>%
    filter(tissue == !!tissue, group == "2+") %>%
    pull(observed)
  
  observed_total <- sum(obs_vals)
  
  sim_totals <- colSums(sim_mat[rownames(sim_mat) == "2+", , drop = FALSE])
  
  p_emp <- mean(sim_totals >= observed_total)
  
  tibble(
    tissue = tissue,
    observed_shared_deg_ge2 = observed_total,
    empirical_p = p_emp
  )
})

pval_table

result <- out$result
plot_df <- tibble(
  Tissue = result$tissue,
  GroupCount = result$group,
  Observed = as.numeric(result$observed),
  Expected = as.numeric(result$expected),
  SD = as.numeric(result$sd)
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD,
    GroupCount_num = ifelse(GroupCount == "1", 1, 2)
  )

pval_positions <- plot_df %>%
  group_by(Tissue) %>%
  summarise(
    # x_pos = max(GroupCount_num, na.rm = TRUE)
    x_pos = 2.5) %>%
  left_join(pval_table, by = c("Tissue" = "tissue")) %>%
  mutate(
    label = paste0("p = ", formatC(empirical_p, format = "e", digits = 2),
                   "\n(≥2 groups,\n1,000 perms)"),
    y_pos = max(plot_df$Expected) * 0.8 # position the label above x=2
  )

plot_long <- plot_df %>%
  dplyr::select(Tissue, GroupCount, Observed, Expected, Lower, Upper) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value") %>%
  mutate(GroupCount_num = ifelse(GroupCount == "1", 1, 2),
         Lower = ifelse(Type == "Observed", NA, Lower),
         Upper = ifelse(Type == "Observed", NA, Upper),
         Type = factor(Type, levels = c("Observed", "Expected")))


ggplot(plot_long, aes(x = GroupCount_num, y = Value, fill = Type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(x = GroupCount_num, ymin = Lower, ymax = Upper, group = Type),
    position = position_dodge(width = 0.7),
    width = 0.2,
    color = "gray20"
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("1", "2+")) +
  scale_fill_manual(values = c("Observed" = "red", "Expected" = "gray80")) +
  facet_grid(~ Tissue, scales = "free_x", space = "free") +
  labs(
    x = "Number of shared cell types",
    y = "Proportion of genes"
  ) +
  geom_text(data = pval_positions,
            aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE,
            size = 3.5, hjust = 1)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.tissue_wise.obs_exp.expr_match.png", width=1100/96, height=220/96, dpi=300)




result$label = ifelse(is.infinite(result$z_score),
                      ifelse(result$z_score > 0, "Inf", "-Inf"),
                      round(result$z_score, 1))
result$y = ifelse(is.infinite(result$z_score),
                  ifelse(result$z_score > 0, max(result$z_score[is.finite(result$z_score)], na.rm = TRUE) + 100,
                         min(result$z_score[is.finite(result$z_score)], na.rm = TRUE) - 10),
                  result$z_score)
ggplot(result, aes(x = factor(n_celltypes_shared), y = z_score)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "lightgrey", fill = "lightgrey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylim(-50, 1400) +
  labs(
    x = "Number of shared cell types",
    y = "Z-score of DEG sharing\n(obs/exp)",
    fill = "Strain"
  ) +
  geom_text(data = result,
            aes(label = label, y = y),
            vjust = ifelse(result$z_score < 0, 1.1, -0.5),
            size = 3) +
  facet_grid(~tissue, scales = "free_x", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.tissue_wise.w_strain.zscore.barplot.expr_match.png", width=1100/96, height=260/96, dpi=300)








analyze_deg_overlap_within_tissue <- function(data,
                                              strain_filter = c("C57BL/6", "SS", "SHR"), 
                                              pval_thresh = 0.05,
                                              logfc_thresh = 0.5) {
  tissue_list <- unique(data$tissue)
  
  result_list <- list()
  
  for (tissue in tissue_list) {
    df <- data %>%
      filter(tissue == !!tissue,
             p_val_adj < pval_thresh,
             abs(avg_log2FC) > logfc_thresh,
             strain %in% strain_filter) %>%
    distinct(gene_name, cell_type)
    
    gene_cell_table <- table(df$gene_name)
    obs_freq <- table(gene_cell_table)
    
    result_df <- tibble(
      tissue = tissue,
      group = names(obs_freq),
      observed = obs_freq
    )
    
    result_list[[tissue]] <- result_df
  }
  
  return(list(
    result = bind_rows(result_list)
  ))
}

out <- analyze_deg_overlap_within_tissue(all_genes)


result <- out$result
plot_df <- tibble(
  Tissue = result$tissue,
  GroupCount = as.factor(as.numeric(result$group)),
  Observed = as.numeric(result$observed)
)
ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray80")) +
  geom_text(aes(y = Observed, label = Observed),
            vjust = -0.5, size = 3, color = "black") +
  labs(
    x = "Number of shared cell types",
    y = "Number of genes"
  ) +
  ylim(0, 1300) +
  facet_grid(~Tissue, scales = "free_x", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.tissue_wise.w_strain.png", width=1000/96, height=250/96, dpi=300)
