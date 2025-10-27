#Load in required packages
library(tidyverse)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)
library(forstringr)
library(ggpubr)
library(ggpmisc)

#Change the path to your directory
path <- "your_path"

#Read in dataframes
cpop_data_WT <- read.csv(file.path(path, "CPOP_data_WT.csv"))
cpop_data_ubr1 <- read.csv(file.path(path, "CPOP_data_ubr1.csv"))
cpop_data_san1 <- read.csv(file.path(path, "CPOP_data_san1.csv"))
cpop_data_original <- read.csv(file.path(path, "CPOP_data_original_2024.csv"))
tile1_count_WT <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/WT_per_tile_variant_counts_tile1.csv"))
tile2_count_WT <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/WT_per_tile_variant_counts_tile2.csv"))
tile3_count_WT <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/WT_per_tile_variant_counts_tile3.csv"))
tile4_count_WT <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/WT_per_tile_variant_counts_tile4.csv"))
tile5_count_WT <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/WT_per_tile_variant_counts_tile5.csv"))
tile1_count_san1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/San1_per_tile_variant_counts_tile1.csv"))
tile2_count_san1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/San1_per_tile_variant_counts_tile2.csv"))
tile3_count_san1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/San1_per_tile_variant_counts_tile3.csv"))
tile4_count_san1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/San1_per_tile_variant_counts_tile4.csv"))
tile5_count_san1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/San1_per_tile_variant_counts_tile5.csv"))
tile1_count_ubr1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/Ubr1_per_tile_variant_counts_tile1.csv"))
tile2_count_ubr1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/Ubr1_per_tile_variant_counts_tile2.csv"))
tile3_count_ubr1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/Ubr1_per_tile_variant_counts_tile3.csv"))
tile4_count_ubr1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/Ubr1_per_tile_variant_counts_tile4.csv"))
tile5_count_ubr1 <- read.csv(file.path(path, "count correlation_WT_san1_ubr1/Ubr1_per_tile_variant_counts_tile5.csv"))


#Fig. S1
#Sequencing count correlation between replicates across tiles and conditions after selection for WT strain
#Tile 1
tile1_rep1_WT <- tile1_count_WT%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile1_rep2_WT <- tile1_count_WT%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile1_rep3_WT <- tile1_count_WT%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile1_merged_WT <- cbind(tile1_rep1_WT, tile1_rep2_WT, tile1_rep3_WT)
colnames(tile1_merged_WT) <- make.unique(names(tile1_merged_WT)) 
tile1_merged_WT_plot <- tile1_merged_WT %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile1_rep12_wt.pdf", ggplot(tile1_merged_WT_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile1_rep13_wt.pdf", ggplot(tile1_merged_WT_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile1_rep23_wt.pdf", ggplot(tile1_merged_WT_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 2
tile2_rep1_WT <- tile2_count_WT%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile2_rep2_WT <- tile2_count_WT%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile2_rep3_WT <- tile2_count_WT%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile2_merged_WT <- cbind(tile2_rep1_WT, tile2_rep2_WT, tile2_rep3_WT)
colnames(tile2_merged_WT) <- make.unique(names(tile2_merged_WT)) 
tile2_merged_WT_plot <- tile2_merged_WT %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile2_rep12_wt.pdf", ggplot(tile2_merged_WT_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile2_rep13_wt.pdf", ggplot(tile2_merged_WT_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile2_rep23_wt.pdf", ggplot(tile2_merged_WT_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 3
tile3_rep1_WT <- tile3_count_WT%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile3_rep2_WT <- tile3_count_WT%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile3_rep3_WT <- tile3_count_WT%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile3_merged_WT <- cbind(tile3_rep1_WT, tile3_rep2_WT, tile3_rep3_WT)
colnames(tile3_merged_WT) <- make.unique(names(tile3_merged_WT)) 
tile3_merged_WT_plot <- tile3_merged_WT %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile3_rep12_wt.pdf", ggplot(tile3_merged_WT_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile3_rep13_wt.pdf", ggplot(tile3_merged_WT_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile3_rep23_wt.pdf", ggplot(tile3_merged_WT_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 4
tile4_rep1_WT <- tile4_count_WT%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile4_rep2_WT <- tile4_count_WT%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile4_rep3_WT <- tile4_count_WT%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile4_merged_WT <- cbind(tile4_rep1_WT, tile4_rep2_WT, tile4_rep3_WT)
colnames(tile4_merged_WT) <- make.unique(names(tile4_merged_WT)) 
tile4_merged_WT_plot <- tile4_merged_WT %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile4_rep12_wt.pdf", ggplot(tile4_merged_WT_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile4_rep13_wt.pdf", ggplot(tile4_merged_WT_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile4_rep23_wt.pdf", ggplot(tile4_merged_WT_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 5
tile5_rep1_WT <- tile5_count_WT%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile5_rep2_WT <- tile5_count_WT%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile5_rep3_WT <- tile5_count_WT%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile5_merged_WT <- cbind(tile5_rep1_WT, tile5_rep2_WT, tile5_rep3_WT)
colnames(tile5_merged_WT) <- make.unique(names(tile5_merged_WT)) 
tile5_merged_WT_plot <- tile5_merged_WT %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile5_rep12_wt.pdf", ggplot(tile5_merged_WT_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile5_rep13_wt.pdf", ggplot(tile5_merged_WT_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile5_rep23_wt.pdf", ggplot(tile5_merged_WT_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Sequencing count correlation between replicates across tiles and conditions after selection for Ubr1 strain
#Tile 1
tile1_rep1_ubr1 <- tile1_count_ubr1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile1_rep2_ubr1 <- tile1_count_ubr1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile1_rep3_ubr1 <- tile1_count_ubr1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile1_merged_ubr1 <- cbind(tile1_rep1_ubr1, tile1_rep2_ubr1, tile1_rep3_ubr1)
colnames(tile1_merged_ubr1) <- make.unique(names(tile1_merged_ubr1)) 
tile1_merged_ubr1_plot <- tile1_merged_ubr1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile1_rep12_ubr1.pdf", ggplot(tile1_merged_ubr1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile1_rep13_ubr1.pdf", ggplot(tile1_merged_ubr1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile1_rep23_ubr1.pdf", ggplot(tile1_merged_ubr1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 2
tile2_rep1_ubr1 <- tile2_count_ubr1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile2_rep2_ubr1 <- tile2_count_ubr1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile2_rep3_ubr1 <- tile2_count_ubr1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30))
tile2_rep3_ubr1_new <- bind_rows(tile2_rep3_ubr1, as_tibble(matrix(NA, nrow = 150, ncol = 0)))
tile2_merged_ubr1 <- cbind(tile2_rep1_ubr1, tile2_rep2_ubr1, tile2_rep3_ubr1_new)
colnames(tile2_merged_ubr1) <- make.unique(names(tile2_merged_ubr1)) 
tile2_merged_ubr1_plot <- tile2_merged_ubr1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35", "35 °C", NA))))
ggsave("tile2_rep12_ubr1.pdf", ggplot(tile2_merged_ubr1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile2_rep13_ubr1.pdf", ggplot(tile2_merged_ubr1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile2_rep23_ubr1.pdf", ggplot(tile2_merged_ubr1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 3
tile3_rep1_ubr1 <- tile3_count_ubr1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile3_rep2_ubr1 <- tile3_count_ubr1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile3_rep3_ubr1 <- tile3_count_ubr1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30))
tile3_rep3_ubr1_new <- bind_rows(tile3_rep3_ubr1, as_tibble(matrix(NA, nrow = 155, ncol = 0)))
tile3_merged_ubr1 <- cbind(tile3_rep1_ubr1, tile3_rep2_ubr1, tile3_rep3_ubr1_new)
colnames(tile3_merged_ubr1) <- make.unique(names(tile3_merged_ubr1)) 
tile3_merged_ubr1_plot <- tile3_merged_ubr1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35", "35 °C", NA))))
ggsave("tile3_rep12_ubr1.pdf", ggplot(tile3_merged_ubr1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile3_rep13_ubr1.pdf", ggplot(tile3_merged_ubr1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile3_rep23_ubr1.pdf", ggplot(tile3_merged_ubr1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 4
tile4_rep1_ubr1 <- tile4_count_ubr1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile4_rep2_ubr1 <- tile4_count_ubr1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile4_rep3_ubr1 <- tile4_count_ubr1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile4_merged_ubr1 <- cbind(tile4_rep1_ubr1, tile4_rep2_ubr1, tile4_rep3_ubr1)
colnames(tile4_merged_ubr1) <- make.unique(names(tile4_merged_ubr1)) 
tile4_merged_ubr1_plot <- tile4_merged_ubr1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile4_rep12_ubr1.pdf", ggplot(tile4_merged_ubr1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile4_rep13_ubr1.pdf", ggplot(tile4_merged_ubr1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile4_rep23_ubr1.pdf", ggplot(tile4_merged_ubr1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 5
tile5_rep1_ubr1 <- tile5_count_ubr1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile5_rep2_ubr1 <- tile5_count_ubr1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile5_rep3_ubr1 <- tile5_count_ubr1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile5_merged_ubr1 <- cbind(tile5_rep1_ubr1, tile5_rep2_ubr1, tile5_rep3_ubr1)
colnames(tile5_merged_ubr1) <- make.unique(names(tile5_merged_ubr1)) 
tile5_merged_ubr1_plot <- tile5_merged_ubr1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile5_rep12_ubr1.pdf", ggplot(tile5_merged_ubr1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile5_rep13_ubr1.pdf", ggplot(tile5_merged_ubr1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile5_rep23_ubr1.pdf", ggplot(tile5_merged_ubr1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Sequencing count correlation between replicates across tiles and conditions after selection for San1 strain
#Tile 1
tile1_rep1_san1 <- tile1_count_san1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile1_rep2_san1 <- tile1_count_san1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile1_rep3_san1 <- tile1_count_san1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile1_merged_san1 <- cbind(tile1_rep1_san1, tile1_rep2_san1, tile1_rep3_san1)
colnames(tile1_merged_san1) <- make.unique(names(tile1_merged_san1)) 
tile1_merged_san1_plot <- tile1_merged_san1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile1_rep12_san1.pdf", ggplot(tile1_merged_san1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile1_rep13_san1.pdf", ggplot(tile1_merged_san1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile1_rep23_san1.pdf", ggplot(tile1_merged_san1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 2
tile2_rep1_san1 <- tile2_count_san1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile2_rep2_san1 <- tile2_count_san1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile2_rep3_san1 <- tile2_count_san1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile2_merged_san1 <- cbind(tile2_rep1_san1, tile2_rep2_san1, tile2_rep3_san1)
colnames(tile2_merged_san1) <- make.unique(names(tile2_merged_san1)) 
tile2_merged_san1_plot <- tile2_merged_san1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile2_rep12_san1.pdf", ggplot(tile2_merged_san1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile2_rep13_san1.pdf", ggplot(tile2_merged_san1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile2_rep23_san1.pdf", ggplot(tile2_merged_san1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 3
tile3_rep1_san1 <- tile3_count_san1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile3_rep2_san1 <- tile3_count_san1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile3_rep3_san1 <- tile3_count_san1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile3_merged_san1 <- cbind(tile3_rep1_san1, tile3_rep2_san1, tile3_rep3_san1)
colnames(tile3_merged_san1) <- make.unique(names(tile3_merged_san1)) 
tile3_merged_san1_plot <- tile3_merged_san1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile3_rep12_san1.pdf", ggplot(tile3_merged_san1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile3_rep13_san1.pdf", ggplot(tile3_merged_san1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile3_rep23_san1.pdf", ggplot(tile3_merged_san1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 4
tile4_rep1_san1 <- tile4_count_san1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile4_rep2_san1 <- tile4_count_san1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile4_rep3_san1 <- tile4_count_san1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile4_merged_san1 <- cbind(tile4_rep1_san1, tile4_rep2_san1, tile4_rep3_san1)
colnames(tile4_merged_san1) <- make.unique(names(tile4_merged_san1)) 
tile4_merged_san1_plot <- tile4_merged_san1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile4_rep12_san1.pdf", ggplot(tile4_merged_san1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile4_rep13_san1.pdf", ggplot(tile4_merged_san1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile4_rep23_san1.pdf", ggplot(tile4_merged_san1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))

#Tile 5
tile5_rep1_san1 <- tile5_count_san1%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no_sele, count_rep1_30, count_rep1_35))
tile5_rep2_san1 <- tile5_count_san1%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no_sele, count_rep2_30, count_rep2_35))
tile5_rep3_san1 <- tile5_count_san1%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no_sele, count_rep3_30, count_rep3_35))
tile5_merged_san1 <- cbind(tile5_rep1_san1, tile5_rep2_san1, tile5_rep3_san1)
colnames(tile5_merged_san1) <- make.unique(names(tile5_merged_san1)) 
tile5_merged_san1_plot <- tile5_merged_san1 %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no_sele") &
                              str_detect(experiment_rep2, "no_sele") &
                              str_detect(experiment_rep3, "no_sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_right(experiment_rep1, 2) == "35" & 
                                            str_right(experiment_rep2, 2) == "35" &
                                            str_right(experiment_rep3, 2) == "35", "35 °C", NA))))
ggsave("tile5_rep12_san1.pdf", ggplot(tile5_merged_san1_plot, aes(count_rep1, count_rep2, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) + 
         labs(x = "Replicate 1", y = "Replicate 2") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile5_rep13_san1.pdf", ggplot(tile5_merged_san1_plot, aes(count_rep1, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 1", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))
ggsave("tile5_rep23_san1.pdf", ggplot(tile5_merged_san1_plot, aes(count_rep2, count_rep3, color = Condition)) +
         geom_point() +
         scale_y_log10() +
         scale_x_log10() +
         stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
         stat_poly_line(show.legend = F) +
         labs(x = "Replicate 2", y = "Replicate 3") + 
         theme_bw() +
         theme(aspect.ratio = 1))


#Fig. S2
#Delta Ubr1/San1 ((Ubr1-San1) score for insertions
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

#Delta Ubr1/San1 ((Ubr1-San1) score for deletions
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