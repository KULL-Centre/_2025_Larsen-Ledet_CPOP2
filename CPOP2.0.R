#Load in required packages
library(tidyverse)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)
library(ggvenn)
library(ggpubr)
library(readxl)

#Change the path to your directory
path <- "your_path"

#Read in dataframes
cpop_data_WT <- read.csv(file.path(path, "CPOP_data_WT.csv"))
cpop_data_ubr1 <- read.csv(file.path(path, "CPOP_data_ubr1.csv"))
cpop_data_san1 <- read.csv(file.path(path, "CPOP_data_san1.csv"))
cpop_data_original <- read.csv(file.path(path, "CPOP_data_original_2024.csv"))
rSASA <- read.csv(file.path(path, "rSASA.csv"), sep = ",")
nanoDSF_rep1 <- read_excel(file.path(path, "nanoDSF_rep1.xlsx")) %>% mutate(Temp = round(Temp, 2), rep = 1)
nanoDSF_rep2 <- read_excel(file.path(path, "nanoDSF_rep2.xlsx")) %>% mutate(Temp = round(Temp, 2), rep = 2)
nanoDSF_rep3 <- read_excel(file.path(path, "nanoDSF_rep3.xlsx")) %>% mutate(Temp = round(Temp, 2), rep = 3)

#Fig. 2
#Heatmap for WT strain for insertions at 30 °C and 35 °C
cpop_ins_WT <- cpop_data_WT %>%
  filter(mut_type == "ins") %>%
  select(pos, score_30, score_35)
cpop_ins_compl_WT <- t(cpop_ins_WT)
cpop_ins_compl_heatmap_WT <- cpop_ins_compl_WT[-1,]
colnames(cpop_ins_compl_heatmap_WT) <- cpop_ins_compl_WT[1,]
WT_heat_ins <- draw(Heatmap(cpop_ins_compl_heatmap_WT,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        colorRamp2(c(-0.36, 1 , 2.36), c("darkred", "white", "darkblue")),
                        name = ("CPOP score"),
                        heatmap_legend_param = list(direction = "horizontal",
                                                    title_position = "topcenter", 
                                                    legend_width = unit(6, "cm"),
                                                    labels_gp = gpar(fontsize = 6),
                                                    title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                    border = "black", 
                                                    at = c(0, 1, 2),
                                                    labels = c("0", "1", "2")),
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_side = "left",
                        border_gp = gpar(col = "black", lwd = 1.5),
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "dimgrey"),
                heatmap_legend_side = "bottom", 
                annotation_legend_side = "bottom",
                annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                     legend_gp = gpar(fill = c("dimgrey")),
                                                     border = c("black", "black"),
                                                     labels_gp = gpar(fontsize = 6),
                                                     row_gap = unit(2, "mm"))))
pdf(file = "WT_heat_ins.pdf", width = 30, height = 6)
WT_heat_ins <- WT_heat_ins
draw(WT_heat_ins)
dev.off()

#Heatmap for ubr1 strain for insertions at 30 °C and 35 °C
cpop_ins_ubr1 <- cpop_data_ubr1 %>%
  filter(mut_type == "ins") %>%
  select(pos, score_30, score_35)
cpop_ins_compl_ubr1 <- t(cpop_ins_ubr1)
cpop_ins_compl_heatmap_ubr1 <- cpop_ins_compl_ubr1[-1,]
colnames(cpop_ins_compl_heatmap_ubr1) <- cpop_ins_compl_ubr1[1,]
ubr1_heat_ins <- draw(Heatmap(cpop_ins_compl_heatmap_ubr1,
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE,
                          colorRamp2(c(-0.30, 1 , 2.30), c("darkred", "white", "darkblue")),
                          name = ("CPOP score"),
                          heatmap_legend_param = list(direction = "horizontal",
                                                      title_position = "topcenter", 
                                                      legend_width = unit(6, "cm"),
                                                      labels_gp = gpar(fontsize = 6),
                                                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                      border = "black", 
                                                      at = c(0, 1, 2),
                                                      labels = c("0", "1", "2")),
                          row_names_gp = gpar(fontsize = 8),
                          column_names_gp = gpar(fontsize = 5),
                          row_names_side = "left",
                          border_gp = gpar(col = "black", lwd = 1.5),
                          rect_gp = gpar(col = "black", lwd = 1),
                          na_col = "dimgrey"),
                  heatmap_legend_side = "bottom", 
                  annotation_legend_side = "bottom",
                  annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                       legend_gp = gpar(fill = c("dimgrey")),
                                                       border = c("black", "black"),
                                                       labels_gp = gpar(fontsize = 6),
                                                       row_gap = unit(2, "mm"))))
pdf(file = "ubr1_heat_ins.pdf", width = 30, height = 6)
ubr1_heat_ins <- ubr1_heat_ins
draw(ubr1_heat_ins)
dev.off()

#Heatmap for san1 strain for insertions at 30 °C and 35 °C
cpop_ins_san1 <- cpop_data_san1 %>%
  filter(mut_type == "ins") %>%
  select(pos, score_30, score_35)
cpop_ins_compl_san1 <- t(cpop_ins_san1)
cpop_ins_compl_heatmap_san1 <- cpop_ins_compl_san1[-1,]
colnames(cpop_ins_compl_heatmap_san1) <- cpop_ins_compl_san1[1,]
san1_heat_ins <- draw(Heatmap(cpop_ins_compl_heatmap_san1,
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE,
                          colorRamp2(c(-0.40, 1 , 2.40), c("darkred", "white", "darkblue")),
                          name = ("CPOP score"),
                          heatmap_legend_param = list(direction = "horizontal",
                                                      title_position = "topcenter", 
                                                      legend_width = unit(6, "cm"),
                                                      labels_gp = gpar(fontsize = 6),
                                                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                      border = "black", 
                                                      at = c(0, 1, 2),
                                                      labels = c("0", "1", "2")),
                          row_names_gp = gpar(fontsize = 8),
                          column_names_gp = gpar(fontsize = 5),
                          row_names_side = "left",
                          border_gp = gpar(col = "black", lwd = 1.5),
                          rect_gp = gpar(col = "black", lwd = 1),
                          na_col = "dimgrey"),
                  heatmap_legend_side = "bottom", 
                  annotation_legend_side = "bottom",
                  annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                       legend_gp = gpar(fill = c("dimgrey")),
                                                       border = c("black", "black"),
                                                       labels_gp = gpar(fontsize = 6),
                                                       row_gap = unit(2, "mm"))))
pdf(file = "san1_heat_ins.pdf", width = 30, height = 6)
san1_heat_ins <- san1_heat_ins
draw(san1_heat_ins)
dev.off()

#Heatmap for WT strain for deletions at 30 °C and 35 °C
cpop_del_WT <- cpop_data_WT %>%
  filter(mut_type == "del") %>%
  select(pos, score_30, score_35)
cpop_del_compl_WT <- t(cpop_del_WT)
cpop_del_compl_heatmap_WT <- cpop_del_compl_WT[-1,]
colnames(cpop_del_compl_heatmap_WT) <- cpop_del_compl_WT[1,]
WT_heat_del <- draw(Heatmap(cpop_del_compl_heatmap_WT,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        colorRamp2(c(-0.31, 1 , 2.31), c("darkred", "white", "darkblue")),
                        name = ("CPOP score"),
                        heatmap_legend_param = list(direction = "horizontal",
                                                    title_position = "topcenter", 
                                                    legend_width = unit(6, "cm"),
                                                    labels_gp = gpar(fontsize = 6),
                                                    title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                    border = "black", 
                                                    at = c(0, 1, 2),
                                                    labels = c("0", "1", "2")),
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_side = "left",
                        border_gp = gpar(col = "black", lwd = 1.5),
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "dimgrey"),
                heatmap_legend_side = "bottom", 
                annotation_legend_side = "bottom",
                annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                     legend_gp = gpar(fill = c("dimgrey")),
                                                     border = c("black", "black"),
                                                     labels_gp = gpar(fontsize = 6),
                                                     row_gap = unit(2, "mm"))))
pdf(file = "WT_heat_del.pdf", width = 30, height = 6)
WT_heat_del <- WT_heat_del
draw(WT_heat_del)
dev.off()

#Heatmap for ubr1 strain for deletions at 30 °C and 35 °C
cpop_del_ubr1 <- cpop_data_ubr1 %>%
  filter(mut_type == "del") %>%
  select(pos, score_30, score_35)
cpop_del_compl_ubr1 <- t(cpop_del_ubr1)
cpop_del_compl_heatmap_ubr1 <- cpop_del_compl_ubr1[-1,]
colnames(cpop_del_compl_heatmap_ubr1) <- cpop_del_compl_ubr1[1,]
ubr1_heat_del <- draw(Heatmap(cpop_del_compl_heatmap_ubr1,
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE,
                            colorRamp2(c(-0.25, 1 , 2.25), c("darkred", "white", "darkblue")),
                            name = ("CPOP score"),
                            heatmap_legend_param = list(direction = "horizontal",
                                                        title_position = "topcenter", 
                                                        legend_width = unit(6, "cm"),
                                                        labels_gp = gpar(fontsize = 6),
                                                        title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                        border = "black", 
                                                        at = c(0, 1, 2),
                                                        labels = c("0", "1", "2")),
                            row_names_gp = gpar(fontsize = 8),
                            column_names_gp = gpar(fontsize = 5),
                            row_names_side = "left",
                            border_gp = gpar(col = "black", lwd = 1.5),
                            rect_gp = gpar(col = "black", lwd = 1),
                            na_col = "dimgrey"),
                    heatmap_legend_side = "bottom", 
                    annotation_legend_side = "bottom",
                    annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                         legend_gp = gpar(fill = c("dimgrey")),
                                                         border = c("black", "black"),
                                                         labels_gp = gpar(fontsize = 6),
                                                         row_gap = unit(2, "mm"))))
pdf(file = "ubr1_heat_del.pdf", width = 30, height = 6)
ubr1_heat_del <- ubr1_heat_del
draw(ubr1_heat_del)
dev.off()

#Heatmap for san1 strain for deletions at 30 °C and 35 °C
cpop_del_san1 <- cpop_data_san1 %>%
  filter(mut_type == "del") %>%
  select(pos, score_30, score_35)
cpop_del_compl_san1 <- t(cpop_del_san1)
cpop_del_compl_heatmap_san1 <- cpop_del_compl_san1[-1,]
colnames(cpop_del_compl_heatmap_san1) <- cpop_del_compl_san1[1,]
san1_heat_del <- draw(Heatmap(cpop_del_compl_heatmap_san1,
                              cluster_rows = FALSE, 
                              cluster_columns = FALSE,
                              colorRamp2(c(-0.36, 1 , 2.36), c("darkred", "white", "darkblue")),
                              name = ("CPOP score"),
                              heatmap_legend_param = list(direction = "horizontal",
                                                          title_position = "topcenter", 
                                                          legend_width = unit(6, "cm"),
                                                          labels_gp = gpar(fontsize = 6),
                                                          title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                          border = "black", 
                                                          at = c(0, 1, 2),
                                                          labels = c("0", "1", "2")),
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 5),
                              row_names_side = "left",
                              border_gp = gpar(col = "black", lwd = 1.5),
                              rect_gp = gpar(col = "black", lwd = 1),
                              na_col = "dimgrey"),
                      heatmap_legend_side = "bottom", 
                      annotation_legend_side = "bottom",
                      annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                           legend_gp = gpar(fill = c("dimgrey")),
                                                           border = c("black", "black"),
                                                           labels_gp = gpar(fontsize = 6),
                                                           row_gap = unit(2, "mm"))))
pdf(file = "san1_heat_del.pdf", width = 30, height = 6)
san1_heat_del <- san1_heat_del
draw(san1_heat_del)
dev.off()


#Fig. 3A
#Calculate mean synonymous and nonsense score across all three datasets (WT, ubr1, san1)
datasets_merged <- rbind(cpop_data_WT, cpop_data_ubr1, cpop_data_san1)
syn_merged <- datasets_merged %>%
  filter(mut_type == "syn") %>%
  summarise(mean_syn = ((mean(score_30, na.rm = T))+(mean(score_35, na.rm = T)))/2)
stop_merged <- datasets_merged %>%
  filter(mut_type == "stop") %>%
  summarise(mean_stop = ((mean(score_30, na.rm = T))+(mean(score_35, na.rm = T)))/2)

#Overlapping score distributions for insertions for WT, ubr1 and san1 strains
cpop_ins_gather_WT <- cpop_data_WT %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_35))
cpop_ins_gather_ubr1 <- cpop_data_ubr1 %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_35))
cpop_ins_gather_san1 <- cpop_data_san1 %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_35))

ggsave("overlapping_score_dist_ins.pdf",
       ggplot() + 
         geom_density(cpop_ins_gather_ubr1, mapping = aes(score, fill = "Ubr1", color = "Ubr1"), alpha = 0.5) +
         geom_density(cpop_ins_gather_san1, mapping = aes(score, fill = "San1", color = "San1"), alpha = 0.5) +
         geom_density(cpop_ins_gather_WT, mapping = aes(score, fill = "WT", color = "WT"), alpha = 0.5) +
         geom_segment(cpop_ins_gather_WT, mapping = aes(x = syn_merged$mean_syn, y = 0, xend = syn_merged$mean_syn, yend = 1.55), color = "black", linetype = "longdash") +
         geom_segment(cpop_ins_gather_WT, mapping = aes(x = stop_merged$mean_stop, y = 0, xend = stop_merged$mean_stop, yend = 1.55), color = "red", linetype = "longdash") +
         annotate("label", x = syn_merged$mean_syn, y = 1.58, label = "Mean syn" , color = "black") +
         annotate("label", x = stop_merged$mean_stop, y = 1.58, label = "Mean stop" , color = "red") +
         facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_35 = "35 °C"))) +
         scale_x_continuous(limits = c(-0.4, 1.2), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         scale_colour_manual(name = NULL, values = c("WT" = "black", "San1" = "black", "Ubr1" = "black"), labels = c("WT" = expression(WT~strain), "San1" = expression(San1*Delta~strain), "Ubr1" = expression(Ubr1*Delta~strain)), breaks = c("WT","San1", "Ubr1")) +
         scale_fill_manual(name = NULL, values = adjustcolor(c("darkorange", "purple", "green4"), alpha.f = 0.5), labels = c("WT" = expression(WT~strain), "San1" = expression(San1*Delta~strain), "Ubr1" = expression(Ubr1*Delta~strain)), breaks = c("WT","San1", "Ubr1")) +
         guides(fill = guide_legend(override.aes = list(alpha = c(0.3, 0.3, 0.3)))) +
         xlab("CPOP score") +
         ylab("Count") +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14, hjust = 0), panel.spacing = unit(1, "lines"), text = element_text(size = 14)))

#Overlapping score distributions for deletions for WT, ubr1 and san1 strains
cpop_del_gather_WT <- cpop_data_WT %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_35))
cpop_del_gather_ubr1 <- cpop_data_ubr1 %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_35))
cpop_del_gather_san1 <- cpop_data_san1 %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_35))

ggsave("overlapping_score_dist_del.pdf",
       ggplot() + 
         geom_density(cpop_del_gather_ubr1, mapping = aes(score, fill = "Ubr1", color = "Ubr1"), alpha = 0.5) +
         geom_density(cpop_del_gather_san1, mapping = aes(score, fill = "San1", color = "San1"), alpha = 0.5) +
         geom_density(cpop_del_gather_WT, mapping = aes(score, fill = "WT", color = "WT"), alpha = 0.5) +
         geom_segment(cpop_del_gather_WT, mapping = aes(x = syn_merged$mean_syn, y = 0, xend = syn_merged$mean_syn, yend = 1.90), color = "black", linetype = "longdash") +
         geom_segment(cpop_del_gather_WT, mapping = aes(x = stop_merged$mean_stop, y = 0, xend = stop_merged$mean_stop, yend = 1.90), color = "red", linetype = "longdash") +
         annotate("label", x = syn_merged$mean_syn, y = 1.91, label = "Mean syn" , color = "black") +
         annotate("label", x = stop_merged$mean_stop, y = 1.91, label = "Mean stop" , color = "red") +
         facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_35 = "35 °C"))) +
         scale_x_continuous(limits = c(-0.4, 1.2), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         scale_colour_manual(name = NULL, values = c("WT" = "black", "San1" = "black", "Ubr1" = "black"), labels = c("WT" = expression(WT~strain), "San1" = expression(San1*Delta~strain), "Ubr1" = expression(Ubr1*Delta~strain)), breaks = c("WT","San1", "Ubr1")) +
         scale_fill_manual(name = NULL, values = adjustcolor(c("darkorange", "purple", "green4"), alpha.f = 0.5), labels = c("WT" = expression(WT~strain), "San1" = expression(San1*Delta~strain), "Ubr1" = expression(Ubr1*Delta~strain)), breaks = c("WT","San1", "Ubr1")) +
         guides(fill = guide_legend(override.aes = list(alpha = c(0.3, 0.3, 0.3)))) +
         xlab("CPOP score") +
         ylab("Count") +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14, hjust = 0), panel.spacing = unit(1, "lines"), text = element_text(size = 14)))


#Fig. 3B
#Heatmap for delta ubr1 (WT vs. ubr1) and delta san1 (WT vs. san1) scores for insertions
datasets_merged_ins <- cbind(cpop_ins_WT, cpop_ins_ubr1, cpop_ins_san1)
colnames(datasets_merged_ins)[1:3] <- paste0(colnames(datasets_merged_ins)[1:3], "_WT")
colnames(datasets_merged_ins)[4:6] <- paste0(colnames(datasets_merged_ins)[4:6], "_ubr1")
colnames(datasets_merged_ins)[7:9] <- paste0(colnames(datasets_merged_ins)[7:9], "_san1")
datasets_merged_ins_delta_WT <- datasets_merged_ins %>%
  mutate(dsan1_30 = score_30_san1 - score_30_WT,
         dsan1_35 = score_35_san1 - score_35_WT,
         dubr1_30 = score_30_ubr1 - score_30_WT,
         dubr1_35 = score_35_ubr1 - score_35_WT) %>%
  select(pos_WT, dsan1_30, dsan1_35, dubr1_30, dubr1_35)
datasets_merged_ins_delta_WT_t <- t(datasets_merged_ins_delta_WT)
cpop_data_delta_WT_heat <- datasets_merged_ins_delta_WT_t[-1, , drop = FALSE]
colnames(cpop_data_delta_WT_heat) <- datasets_merged_ins_delta_WT_t[1,]
delta_WT_ins <- draw(Heatmap(cpop_data_delta_WT_heat,
                                 cluster_rows = FALSE, 
                                 cluster_columns = FALSE,
                                 colorRamp2(c(0, 0.48), c("white", "magenta")),
                                 heatmap_legend_param = list(direction = "horizontal",
                                                             title = expression(Delta*"CPOP score"),
                                                             title_position = "topcenter", 
                                                             legend_width = unit(6, "cm"),
                                                             labels_gp = gpar(fontsize = 6),
                                                             title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                             border = "black", 
                                                             at = c(0, 0.5),
                                                             labels = c(0, 0.5)),
                                 row_names_gp = gpar(fontsize = 8),
                                 column_names_gp = gpar(fontsize = 5),
                                 row_names_side = "left",
                                 border_gp = gpar(col = "black", lwd = 1),
                                 rect_gp = gpar(col = "black", lwd = 1),
                                 na_col = "dimgrey"),
                         heatmap_legend_side = "bottom", 
                         annotation_legend_side = "bottom",
                         annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                              legend_gp = gpar(fill = c("dimgrey")),
                                                              border = c("black", "black"),
                                                              labels_gp = gpar(fontsize = 6),
                                                              row_gap = unit(2, "mm"))))
pdf(file = "delta_WT_ins.pdf", width = 30, height = 6)
delta_WT_ins <- delta_WT_ins
draw(delta_WT_ins)
dev.off()

#Heatmap for delta ubr1 (WT vs. ubr1) and delta san1 (WT vs. san1) scores for deletions
datasets_merged_del <- cbind(cpop_del_WT, cpop_del_ubr1, cpop_del_san1)
colnames(datasets_merged_del)[1:3] <- paste0(colnames(datasets_merged_del)[1:3], "_WT")
colnames(datasets_merged_del)[4:6] <- paste0(colnames(datasets_merged_del)[4:6], "_ubr1")
colnames(datasets_merged_del)[7:9] <- paste0(colnames(datasets_merged_del)[7:9], "_san1")
datasets_merged_del_delta_WT <- datasets_merged_del %>%
  mutate(dsan1_30 = score_30_san1 - score_30_WT,
         dsan1_35 = score_35_san1 - score_35_WT,
         dubr1_30 = score_30_ubr1 - score_30_WT,
         dubr1_35 = score_35_ubr1 - score_35_WT) %>%
  select(pos_WT, dsan1_30, dsan1_35, dubr1_30, dubr1_35)
datasets_merged_del_delta_WT_t <- t(datasets_merged_del_delta_WT)
cpop_data_delta_WT_heat <- datasets_merged_del_delta_WT_t[-1, , drop = FALSE]
colnames(cpop_data_delta_WT_heat) <- datasets_merged_del_delta_WT_t[1,]
delta_WT_del <- draw(Heatmap(cpop_data_delta_WT_heat,
                             cluster_rows = FALSE, 
                             cluster_columns = FALSE,
                             colorRamp2(c(0, 0.53), c("white", "magenta")),
                             heatmap_legend_param = list(direction = "horizontal",
                                                         title = expression(Delta*"CPOP score"),
                                                         title_position = "topcenter", 
                                                         legend_width = unit(6, "cm"),
                                                         labels_gp = gpar(fontsize = 6),
                                                         title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                         border = "black", 
                                                         at = c(0, 0.5),
                                                         labels = c(0, 0.5)),
                             row_names_gp = gpar(fontsize = 8),
                             column_names_gp = gpar(fontsize = 5),
                             row_names_side = "left",
                             border_gp = gpar(col = "black", lwd = 1),
                             rect_gp = gpar(col = "black", lwd = 1),
                             na_col = "dimgrey"),
                     heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom",
                     annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                          legend_gp = gpar(fill = c("dimgrey")),
                                                          border = c("black", "black"),
                                                          labels_gp = gpar(fontsize = 6),
                                                          row_gap = unit(2, "mm"))))
pdf(file = "delta_WT_del.pdf", width = 30, height = 6)
delta_WT_del <- delta_WT_del
draw(delta_WT_del)
dev.off()


#Fig. 4A
#Group secondary structure elements
rSASA_group <- rSASA %>%
  mutate(pos = DSSP.index + 1,
         sec = ifelse(Secondary.structure == "G" | 
                        Secondary.structure == "H" | 
                        Secondary.structure == "I", "helices",
                      ifelse(Secondary.structure == "E" | 
                               Secondary.structure == "B", "sheets",
                             ifelse(Secondary.structure == "T" |
                                      Secondary.structure == "S" |
                                      Secondary.structure == "-", "loops", NA))))

#Box plots for secondary structure elements for WT strain for insertions
rSASA_cpop_WT_ins <- merge(cpop_data_WT, rSASA_group, by = "pos") %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_35))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
rSASA_cpop_WT_ins$sec <- factor(rSASA_cpop_WT_ins$sec, levels = c("loops", "helices", "sheets"))
ggsave("rSASA_cpop_WT_ins.pdf", ggplot(rSASA_cpop_WT_ins, aes(experiment, score, fill = sec), na.rm = F) +
         stat_boxplot(geom ='errorbar') +
         geom_boxplot() +
         scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                           labels=c("Loops", "Helices", "Sheets"),
                           name = "Secondary structure") +
         scale_x_discrete(labels = c("30 °C", "35 °C")) +
         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0))

#Box plots for secondary structure elements for WT strain for deletions
rSASA_cpop_WT_del <- merge(cpop_data_WT, rSASA_group, by = "pos") %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_35))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
rSASA_cpop_WT_del$sec <- factor(rSASA_cpop_WT_del$sec, levels = c("loops", "helices", "sheets"))
ggsave("rSASA_cpop_WT_del.pdf", ggplot(rSASA_cpop_WT_del, aes(experiment, score, fill = sec), na.rm = F) +
         stat_boxplot(geom ='errorbar') +
         geom_boxplot() +
         scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                           labels=c("Loops", "Helices", "Sheets"),
                           name = "Secondary structure") +
         scale_x_discrete(labels = c("30 °C", "35 °C")) +
         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0))

#Box plots for secondary structure elements for ubr1 strain for insertions
rSASA_cpop_ubr1_ins <- merge(cpop_data_ubr1, rSASA_group, by = "pos") %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_35))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
rSASA_cpop_ubr1_ins$sec <- factor(rSASA_cpop_ubr1_ins$sec, levels = c("loops", "helices", "sheets"))
ggsave("rSASA_cpop_ubr1_ins.pdf", ggplot(rSASA_cpop_ubr1_ins, aes(experiment, score, fill = sec), na.rm = F) +
         stat_boxplot(geom ='errorbar') +
         geom_boxplot() +
         scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                           labels=c("Loops", "Helices", "Sheets"),
                           name = "Secondary structure") +
         scale_x_discrete(labels = c("30 °C", "35 °C")) +
         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0))

#Box plots for secondary structure elements for ubr1 strain for deletions
rSASA_cpop_ubr1_del <- merge(cpop_data_ubr1, rSASA_group, by = "pos") %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_35))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
rSASA_cpop_ubr1_del$sec <- factor(rSASA_cpop_ubr1_del$sec, levels = c("loops", "helices", "sheets"))
ggsave("rSASA_cpop_ubr1_del.pdf", ggplot(rSASA_cpop_ubr1_del, aes(experiment, score, fill = sec), na.rm = F) +
         stat_boxplot(geom ='errorbar') +
         geom_boxplot() +
         scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                           labels=c("Loops", "Helices", "Sheets"),
                           name = "Secondary structure") +
         scale_x_discrete(labels = c("30 °C", "35 °C")) +
         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0))

#Box plots for secondary structure elements for san1 strain for insertions
rSASA_cpop_san1_ins <- merge(cpop_data_san1, rSASA_group, by = "pos") %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_35))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
rSASA_cpop_san1_ins$sec <- factor(rSASA_cpop_san1_ins$sec, levels = c("loops", "helices", "sheets"))
ggsave("rSASA_cpop_san1_ins.pdf", ggplot(rSASA_cpop_san1_ins, aes(experiment, score, fill = sec), na.rm = F) +
         stat_boxplot(geom ='errorbar') +
         geom_boxplot() +
         scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                           labels=c("Loops", "Helices", "Sheets"),
                           name = "Secondary structure") +
         scale_x_discrete(labels = c("30 °C", "35 °C")) +
         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0))

#Box plots for secondary structure elements for san1 strain for deletions
rSASA_cpop_san1_del <- merge(cpop_data_san1, rSASA_group, by = "pos") %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_35))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
rSASA_cpop_san1_del$sec <- factor(rSASA_cpop_san1_del$sec, levels = c("loops", "helices", "sheets"))
ggsave("rSASA_cpop_san1_del.pdf", ggplot(rSASA_cpop_san1_del, aes(experiment, score, fill = sec), na.rm = F) +
         stat_boxplot(geom ='errorbar') +
         geom_boxplot() +
         scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                           labels=c("Loops", "Helices", "Sheets"),
                           name = "Secondary structure") +
         scale_x_discrete(labels = c("30 °C", "35 °C")) +
         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0))


#Fig. 4B
#Delta Ubr1/San1 (Ubr1-San1) score for insertions
cpop_ins_WT <- cpop_data_WT %>%
  filter(mut_type == "ins") %>%
  select(pos, score_30, score_35)
cpop_ins_ubr1 <- cpop_data_ubr1 %>%
  filter(mut_type == "ins") %>%
  select(pos, score_30, score_35)
cpop_ins_san1 <- cpop_data_san1 %>%
  filter(mut_type == "ins") %>%
  select(pos, score_30, score_35)
datasets_merged_ins <- cbind(cpop_ins_WT, cpop_ins_ubr1, cpop_ins_san1)
colnames(datasets_merged_ins)[1:3] <- paste0(colnames(datasets_merged_ins)[1:3], "_WT")
colnames(datasets_merged_ins)[4:6] <- paste0(colnames(datasets_merged_ins)[4:6], "_ubr1")
colnames(datasets_merged_ins)[7:9] <- paste0(colnames(datasets_merged_ins)[7:9], "_san1")
datasets_merged_ins_delta_strain <- datasets_merged_ins %>%
  mutate("dstrain_30" = score_30_ubr1 - score_30_san1,
         "dstrain_35" = score_35_ubr1 - score_35_san1) %>%
  select(pos_WT, "dstrain_30", "dstrain_35")
datasets_merged_ins_delta_strain_t <- t(datasets_merged_ins_delta_strain)
cpop_data_delta_strain_heat <- datasets_merged_ins_delta_strain_t[-1, , drop = FALSE]
colnames(cpop_data_delta_strain_heat) <- datasets_merged_ins_delta_strain_t[1,]
delta_strain_ins <- draw(Heatmap(cpop_data_delta_strain_heat,
                                 cluster_rows = FALSE, 
                                 cluster_columns = FALSE,
                                 colorRamp2(c(-0.25, 0, 0.25), c("cyan", "white", "magenta")),
                                 heatmap_legend_param = list(direction = "horizontal",
                                                             title = expression(Delta*"strain score"),
                                                             title_position = "topcenter", 
                                                             legend_width = unit(6, "cm"),
                                                             labels_gp = gpar(fontsize = 6),
                                                             title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                             border = "black", 
                                                             at = c(-0.25, 0, 0.25),
                                                             labels = c("-0.25\nSan1\nspecific", "0\nEqually\nspecific", "0.25\nUbr1\nspecific")),
                                 row_names_gp = gpar(fontsize = 8),
                                 column_names_gp = gpar(fontsize = 5),
                                 row_names_side = "left",
                                 border_gp = gpar(col = "black", lwd = 1),
                                 rect_gp = gpar(col = "black", lwd = 1),
                                 na_col = "dimgrey"),
                         heatmap_legend_side = "bottom", 
                         annotation_legend_side = "bottom",
                         annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                              legend_gp = gpar(fill = c("dimgrey")),
                                                              border = c("black", "black"),
                                                              labels_gp = gpar(fontsize = 6),
                                                              row_gap = unit(2, "mm"))))
pdf(file = "delta_strain_ins.pdf", width = 30, height = 6)
delta_strain_ins <- delta_strain_ins
draw(delta_strain_ins)
dev.off()

#Delta Ubr1/San1 (Ubr1-San1) score for deletions
cpop_del_WT <- cpop_data_WT %>%
  filter(mut_type == "del") %>%
  select(pos, score_30, score_35)
cpop_del_ubr1 <- cpop_data_ubr1 %>%
  filter(mut_type == "del") %>%
  select(pos, score_30, score_35)
cpop_del_san1 <- cpop_data_san1 %>%
  filter(mut_type == "del") %>%
  select(pos, score_30, score_35)
datasets_merged_del <- cbind(cpop_del_WT, cpop_del_ubr1, cpop_del_san1)
colnames(datasets_merged_del)[1:3] <- paste0(colnames(datasets_merged_del)[1:3], "_WT")
colnames(datasets_merged_del)[4:6] <- paste0(colnames(datasets_merged_del)[4:6], "_ubr1")
colnames(datasets_merged_del)[7:9] <- paste0(colnames(datasets_merged_del)[7:9], "_san1")
datasets_merged_del_delta_strain <- datasets_merged_del %>%
  mutate("dstrain_30" = score_30_ubr1 - score_30_san1,
         "dstrain_35" = score_35_ubr1 - score_35_san1) %>%
  select(pos_WT, "dstrain_30", "dstrain_35")
datasets_merged_del_delta_strain_t <- t(datasets_merged_del_delta_strain)
cpop_data_delta_strain_heat <- datasets_merged_del_delta_strain_t[-1, , drop = FALSE]
colnames(cpop_data_delta_strain_heat) <- datasets_merged_del_delta_strain_t[1,]
delta_strain_del <- draw(Heatmap(cpop_data_delta_strain_heat,
                                 cluster_rows = FALSE, 
                                 cluster_columns = FALSE,
                                 colorRamp2(c(-0.25, 0, 0.25), c("cyan", "white", "magenta")),
                                 heatmap_legend_param = list(direction = "horizontal",
                                                             title = expression(Delta*"strain score"),
                                                             title_position = "topcenter", 
                                                             legend_width = unit(6, "cm"),
                                                             labels_gp = gpar(fontsize = 6),
                                                             title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                             border = "black", 
                                                             at = c(-0.25, 0, 0.25),
                                                             labels = c("-0.25\nSan1\nspecific", "0\nEqually\nspecific", "0.25\nUbr1\nspecific")),
                                 row_names_gp = gpar(fontsize = 8),
                                 column_names_gp = gpar(fontsize = 5),
                                 row_names_side = "left",
                                 border_gp = gpar(col = "black", lwd = 1),
                                 rect_gp = gpar(col = "black", lwd = 1),
                                 na_col = "dimgrey"),
                         heatmap_legend_side = "bottom", 
                         annotation_legend_side = "bottom",
                         annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                                              legend_gp = gpar(fill = c("dimgrey")),
                                                              border = c("black", "black"),
                                                              labels_gp = gpar(fontsize = 6),
                                                              row_gap = unit(2, "mm"))))
pdf(file = "delta_strain_del.pdf", width = 30, height = 6)
delta_strain_del <- delta_strain_del
draw(delta_strain_del)
dev.off()


#Fig. 7E
#Melting curves
nanoDSF_merged <- bind_rows(nanoDSF_rep1, nanoDSF_rep2, nanoDSF_rep3)
nanoDSF_merged_long <- nanoDSF_merged %>%
  pivot_longer(cols = -c(Temp, rep),
               names_to = "variable",
               values_to = "value") %>%
  mutate(variable = sub("_[123]$", "", variable))
nanoDSF_merged_mean <- nanoDSF_merged_long %>%
  group_by(Temp, variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
nanoDSF_merged_mean_plot <- nanoDSF_merged_mean %>%
  filter(Temp <= 65) %>%
  separate(variable, into = c("variant", "substrate"), sep = "_") %>%
  mutate(mean_value = pmin(mean_value, 1),
         variant = factor(variant, levels = c("WT", "D22del", "E173ins"),
                          labels = c("WT", "D22∆", "E173ins")))
ggsave("Melting_curves.pdf", ggplot(nanoDSF_merged_mean_plot, aes(x = Temp, y = mean_value, color = substrate)) +
  geom_line(linewidth = 2) +
  facet_wrap(~ variant, ncol = 1) +
  scale_color_manual(values = c(folate = "#FFCC66", MTX = "black")) +
  xlab("Temperature (°C)") +
  ylab("Fraction unfolded") +
  theme_bw(),
  device = cairo_pdf, width = 7, height = 5)


#Fig. 8A
#Venn diagram for MTX, TS, ubr1 and san1 for insertions
datasets_merged_ins_delta_WT1 <- datasets_merged_ins_delta_WT %>%
  rename(pos = pos_WT)
cpop_delta_WT_mtx_ubr1_san1_ins <- cpop_data_original %>%
  filter(mut_type == "ins") %>%
  mutate(dmtx_37 = score_37_MTX - score_37) %>%
  mutate(dtemp = score_30 - score_37) %>%
  merge(datasets_merged_ins_delta_WT1, by = "pos") %>%
  select(pos, dtemp, dmtx_37, dubr1_30, dsan1_30)
#Filter top 25% hits by reducing temp, adding MTX, removing ubr1 or removing san1
cpop_delta_WT_mtx_ubr1_san1_ins_25pct_MTX <- cpop_delta_WT_mtx_ubr1_san1_ins %>%
  filter(dmtx_37 >= quantile(dmtx_37, 0.75, na.rm = TRUE))
cpop_delta_WT_mtx_ubr1_san1_ins_25pct_temp <- cpop_delta_WT_mtx_ubr1_san1_ins %>%
  filter(dtemp >= quantile(dtemp, 0.75, na.rm = TRUE))
cpop_delta_WT_mtx_ubr1_san1_ins_25pct_ubr1 <- cpop_delta_WT_mtx_ubr1_san1_ins %>%
  filter(dubr1_30 >= quantile(dubr1_30, 0.75, na.rm = TRUE))
cpop_delta_WT_mtx_ubr1_san1_ins_25pct_san1 <- cpop_delta_WT_mtx_ubr1_san1_ins %>%
  filter(dsan1_30 >= quantile(dsan1_30, 0.75, na.rm = TRUE))
ggsave("venn_ins_25pct.pdf", ggvenn(list(MTX = cpop_delta_WT_mtx_ubr1_san1_ins_25pct_MTX$pos, 
                                         TS = cpop_delta_WT_mtx_ubr1_san1_ins_25pct_temp$pos, 
                                         Ubr1 = cpop_delta_WT_mtx_ubr1_san1_ins_25pct_ubr1$pos, 
                                         San1 = cpop_delta_WT_mtx_ubr1_san1_ins_25pct_san1$pos),
                                    fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")))

#Venn diagram for MTX, TS, ubr1 and san1 for deletions
datasets_merged_del_delta_WT1 <- datasets_merged_del_delta_WT %>%
  rename(pos = pos_WT)
cpop_delta_WT_mtx_ubr1_san1_del <- cpop_data_original %>%
  filter(mut_type == "del") %>%
  mutate(dmtx_37 = score_37_MTX - score_37) %>%
  mutate(dtemp = score_30 - score_37) %>%
  merge(datasets_merged_del_delta_WT1, by = "pos") %>%
  select(pos, dtemp, dmtx_37, dubr1_30, dsan1_30)
#Filter top 25% hits by reducing temp, adding MTX, removing ubr1 or removing san1
cpop_delta_WT_mtx_ubr1_san1_del_25pct_MTX <- cpop_delta_WT_mtx_ubr1_san1_del %>%
  filter(dmtx_37 >= quantile(dmtx_37, 0.75, na.rm = TRUE))
cpop_delta_WT_mtx_ubr1_san1_del_25pct_temp <- cpop_delta_WT_mtx_ubr1_san1_del %>%
  filter(dtemp >= quantile(dtemp, 0.75, na.rm = TRUE))
cpop_delta_WT_mtx_ubr1_san1_del_25pct_ubr1 <- cpop_delta_WT_mtx_ubr1_san1_del %>%
  filter(dubr1_30 >= quantile(dubr1_30, 0.75, na.rm = TRUE))
cpop_delta_WT_mtx_ubr1_san1_del_25pct_san1 <- cpop_delta_WT_mtx_ubr1_san1_del %>%
  filter(dsan1_30 >= quantile(dsan1_30, 0.75, na.rm = TRUE))
ggsave("venn_del_25pct.pdf", ggvenn(list(MTX = cpop_delta_WT_mtx_ubr1_san1_del_25pct_MTX$pos, 
                                         TS = cpop_delta_WT_mtx_ubr1_san1_del_25pct_temp$pos, 
                                         Ubr1 = cpop_delta_WT_mtx_ubr1_san1_del_25pct_ubr1$pos, 
                                         San1 = cpop_delta_WT_mtx_ubr1_san1_del_25pct_san1$pos),
                                    fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")))
