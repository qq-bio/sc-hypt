library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(openxlsx)

rat_df <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Copy of rs064_all_qPCR_data_single_cell_manuscript.xlsx", "rat_mod")
ivsmc_df <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Copy of rs064_all_qPCR_data_single_cell_manuscript.xlsx", "iVSMC_mod")
iec_df <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Copy of rs064_all_qPCR_data_single_cell_manuscript.xlsx", "iEC_mod")
ivsmc_chip_df <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Copy of rs064_all_qPCR_data_single_cell_manuscript.xlsx", "iVC_chIP_mod")
iec_chip_df <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Copy of rs064_all_qPCR_data_single_cell_manuscript.xlsx", "iEC_ChIP_mod")





rat_df_melt <- melt(rat_df, id.vars = "Type", variable.name = "Gene", value.name = "Expression") 
rat_df_melt$Gene <- factor(str_to_title(rat_df_melt$Gene), levels = c("Mrps6", "Slc5a3", "Smim11", "Rcan1"))
d
stat.test <- rat_df_melt %>%
  group_by(Gene) %>%
  t_test(Expression ~ Type) %>%
  add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns" 
  )) %>% 
  filter(p.signif != "ns")

bp <- ggbarplot(
  rat_df_melt, x = "Gene", y = "Expression", add = "mean_sd", 
  fill = "Type",  color = "black", 
  position = position_dodge(0.8)
) +
  scale_fill_manual(values = c("white", "white")) +
  geom_jitter(aes(shape = Type), size = 1.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_shape_manual(values = c(16, 1)) + 
  theme(legend.position = "right") +
  guides(fill = "none") +  
  labs(y = expression("mRNA (copies/" ~ 10^6 ~ " 18S rRNA)"), x = "", shape="")

bp + stat_pvalue_manual(
  stat.test, label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0.5 
)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig5f.png", width=460/96, height=279/96, dpi=300)


# bp <- ggbarplot(
#   rat_df_melt, x = "Gene", y = "Expression", add = "mean_sd", 
#   color= "Type", # palette = c("#00AFBB", "#E7B800"),
#   position = position_dodge(0.8)
# ) +
#   geom_jitter(aes(color = Type), size = 1.5,
#               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
#   theme(legend.position = "right") +
#   labs(y = expression("mRNA (copies/" ~ 10^6 ~ " 18S rRNA)"), x = "", color="")
# 
# stat.test <- rat_df_melt %>%
#   group_by(Gene) %>%
#   t_test(Expression ~ Type) %>%
#   add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
#   mutate(p.signif = case_when(
#     p <= 0.001 ~ "***",
#     p <= 0.01 ~ "**",
#     p <= 0.05 ~ "*",
#     TRUE ~ "ns" 
#   ))
# 
# stat.test <- stat.test %>% filter(p.signif != "ns")
# bp + stat_pvalue_manual(
#   stat.test, label = "p.signif", tip.length = 0.01,
#   bracket.nudge.y = 0.5 
# )
# ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig5f.png", width=460/96, height=279/96, dpi=300)




ivsmc_df_melt <- melt(ivsmc_df, id.vars = "Type", variable.name = "Gene", value.name = "Expression") 
bp <- ggbarplot(
  ivsmc_df_melt, x = "Gene", y = "Expression", add = "mean_sd", 
  color= "Type", palette = c("#D55E00", "#56B4E9"),
  position = position_dodge(0.8)
) +
  geom_jitter(aes(color = Type), size = 1.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(legend.position = "right") +
  labs(y = "Normalized with larger PP allele", x = "", color="")

stat.test <- ivsmc_df_melt %>%
  group_by(Gene) %>%
  t_test(Expression ~ Type) %>%
  add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns" 
  ))

stat.test <- stat.test %>% filter(p.signif != "ns")
bp + stat_pvalue_manual(
  stat.test, label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0.5 
)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig5g.png", width=641/96, height=279/96, dpi=300)





iec_df_melt <- melt(iec_df, id.vars = "Type", variable.name = "Gene", value.name = "Expression") 
bp <- ggbarplot(
  iec_df_melt, x = "Gene", y = "Expression", add = "mean_sd", 
  color= "Type", palette = c("#D55E00", "#56B4E9"),
  position = position_dodge(0.8)
) +
  geom_jitter(aes(color = Type), size = 1.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(legend.position = "right") +
  labs(y = "Normalized with larger PP allele", x = "", color="")

stat.test <- iec_df_melt %>%
  group_by(Gene) %>%
  t_test(Expression ~ Type) %>%
  add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns" 
  ))

stat.test <- stat.test %>% filter(p.signif != "ns")
bp + stat_pvalue_manual(
  stat.test, label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0.5 
)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig5g.iec.png", width=641/96, height=279/96, dpi=300)





H3K4me1_df <- ivsmc_chip_df[, c("H3K4me1", "Type")]
bp <- ggbarplot(
  H3K4me1_df, x = "Type", y = "H3K4me1", add = "mean_sd", 
  color= "Type", palette = c("#D55E00", "#56B4E9"),
  position = position_dodge(0.8)
) +
  geom_jitter(aes(color = Type), size = 1.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  labs(y = "% of input", x = "H3K4me1", color="")

stat.test <- H3K4me1_df %>%
  t_test(H3K4me1 ~ Type, paired = T) %>%
  add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
  mutate(
    xmin = 1, 
    xmax = 2,
    p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns" 
  ))

stat.test <- stat.test %>% filter(p.signif != "ns")
p1 = bp + stat_pvalue_manual(
  stat.test, label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0, inherit.aes = F
)


H3K27ac_df <- ivsmc_chip_df[, c("H3K27ac", "Type")]
bp <- ggbarplot(
  H3K27ac_df, x = "Type", y = "H3K27ac", add = "mean_sd", 
  color= "Type", palette = c("#D55E00", "#56B4E9"),
  position = position_dodge(0.8)
) +
  geom_jitter(aes(color = Type), size = 1.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(legend.position = "right",
        axis.text.x = element_blank()) +
  labs(y = "% of input", x = "H3K27ac", color="")

stat.test <- H3K27ac_df %>%
  t_test(H3K27ac ~ Type, paired = T) %>%
  add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
  mutate(
    xmin = 1, 
    xmax = 2,
    p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns" 
  ))

stat.test <- stat.test %>% filter(p.signif != "ns")
p2 = bp + stat_pvalue_manual(
  stat.test, label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0.5 
)

p1+p2
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig5i.png", width=434/96, height=279/96, dpi=300)








H3K4me1_df <- iec_chip_df[, c("H3K4me1", "Type")]
bp <- ggbarplot(
  H3K4me1_df, x = "Type", y = "H3K4me1", add = "mean_sd", 
  color= "Type", palette = c("#D55E00", "#56B4E9"),
  position = position_dodge(0.8)
) +
  geom_jitter(aes(color = Type), size = 1.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  labs(y = "% of input", x = "H3K4me1", color="")

stat.test <- H3K4me1_df %>%
  t_test(H3K4me1 ~ Type, paired = T) %>%
  add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
  mutate(
    xmin = 1, 
    xmax = 2,
    p.signif = case_when(
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      TRUE ~ "ns" 
    ))

stat.test <- stat.test %>% filter(p.signif != "ns")
p1 = bp + stat_pvalue_manual(
  stat.test, label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0, inherit.aes = F
)


H3K27ac_df <- iec_chip_df[, c("H3K27ac", "Type")]
bp <- ggbarplot(
  H3K27ac_df, x = "Type", y = "H3K27ac", add = "mean_sd", 
  color= "Type", palette = c("#D55E00", "#56B4E9"),
  position = position_dodge(0.8)
) +
  geom_jitter(aes(color = Type), size = 1.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(legend.position = "right",
        axis.text.x = element_blank()) +
  labs(y = "% of input", x = "H3K27ac", color="")

stat.test <- H3K27ac_df %>%
  t_test(H3K27ac ~ Type, paired = T) %>%
  add_xy_position(fun = "mean_sd", x = "Gene", dodge = 0.8) %>%
  mutate(
    xmin = 1, 
    xmax = 2,
    p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns" 
  ))

stat.test <- stat.test %>% filter(p.signif != "ns")
p2 = bp + stat_pvalue_manual(
  stat.test, label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0 
)

p1+p2
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig5i.iEC.png", width=434/96, height=279/96, dpi=300)


