# Figure 5: Single-read transcript classification

library(here)
source(here("scripts/helpers.R"))

# load mapped data ----
total_frame_anno_g <- vroom(here("data/Rdata/total_frame_anno_g.tsv"))

dir <- paste0(here("data"), "/")

## start sites ====
tss <- vroom(paste0(dir, "Rdata/tss_peaks_all_longread.tsv")) %>%
  dplyr::filter(biotype == "mRNA ORF", set == "twt", TSS_class == "primary") %>%
  dplyr::rename(gene = gene_id) %>%
  dplyr::select(gene, UTR5, TSS, cov_peak)

## termination sites ====
tts <- vroom(paste0(dir, "Rdata/tts_peaks_all_longread_cov5.tsv")) %>%
  dplyr::filter(biotype == "mRNA ORF", set == "twt") %>%
  dplyr::rename(gene = gene_id) %>%
  dplyr::select(gene, UTR3, TTS, cov_peak)

total_frame_sites <- total_frame_anno_g %>%
  left_join(tss, by = c("mapped_gene.y" = "gene")) %>%
  left_join(tts, by = c("mapped_gene.y" = "gene"))

total_frame_sites_f <- total_frame_sites %>%
  dplyr::filter(barcode.x %in% c(paste0("barcode",20:23))) %>%
  dplyr::filter(!is.na(TSS), !is.na(TTS)) %>%
  group_by(mapped_gene.y) %>%
  dplyr::mutate(total = n()) %>%
  arrange(desc(total))


total_frame_sites_f2 <- total_frame_sites_f %>%
  rowwise() %>%
  dplyr::filter(start.x %in% TSS:TTS, end.x %in% TSS:TTS) %>%
  ungroup() %>%
  arrange(desc(end.x)) %>%
  dplyr::mutate(rown = 1:nrow(.)) %>%
  dplyr::mutate(category = ifelse(strand.x == "-" & start.x >= (TTS+10) & end.x > (TSS-10), "3_degraded", 
                                  ifelse(strand.x == "-" & end.x <= (TSS-10) & start.x < (TTS+10), "5_degraded",
                                         ifelse(strand.x == "-" & end.x <= (TSS-10) & start.x >= (TTS+10), "degraded",
                                                ifelse(strand.x == "+" & start.x >= (TSS+10) & end.x > (TTS-10), "5_degraded", 
                                                       ifelse(strand.x == "+" & end.x <= (TTS-10) & start.x < (TSS+10), "3_degraded",
                                                              ifelse(strand.x == "+" & end.x <= (TTS-10) & start.x >= (TSS+10), "degraded","full-length"))))))) 

ont_run_info <- read_xlsx(here("data/samples/ont_run_info.xlsx"))

total_frame_sites_f2_quant <- total_frame_sites_f2 %>%
  group_by(mapped_gene.y, condition) %>%
  dplyr::mutate(total = n()) %>%
  group_by(mapped_gene.y, condition,category) %>%
  reframe(frac = n()/total * 100) %>%
  distinct(mapped_gene.y, frac, condition, category)


dir.create(here("figures"), showWarnings = FALSE)
pdf(here("figures/Fig5C_transcript_classes.pdf"), height = 5,width = 5)
p <- ggplot(data = total_frame_sites_f2_quant %>%
         dplyr::mutate(condition = factor(condition, levels = rev(c("WT", "dRNaseR")))) %>%
         dplyr::mutate(category = factor(category, levels = rev(c("full-length", 
                                                              "5_degraded", 
                                                              "3_degraded", 
                                                              "degraded")))),
       aes(x = category, fill = condition, y = frac)) +
  #geom_split_violin(color = NA, alpha = 0.5, scale = "width") +
  geom_jitter(aes(color = condition), position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.6)) +
  geom_boxplot(alpha = 0.6,notchwidth = 0.6,notch = T, color = "black", outlier.colour = NA)+
  theme(axis.text.x = element_text(angle = 90)) +
  theme_pubclean() +
  scale_fill_manual(values = c("#D55E00","#A9A9A9")) +
  scale_color_manual(values = c("#D55E00","#A9A9A9")) +
  ggpubr::stat_compare_means(method = "wilcox", 
                             method.args = list(var.equal = FALSE),
                             label = "p.signif") +
  coord_flip()
print(p)
dev.off()

total_frame_sites_f2_plot <- total_frame_sites_f2 %>%
  ungroup() %>%
  dplyr::select(start.x, end.x, barcode.x, mapped_gene.y, category) %>%
  dplyr::filter(mapped_gene.y %in% "gene438") %>%
  left_join(ont_run_info, by = c("barcode.x" = "barcode")) %>%
  group_by(condition, category, mapped_gene.y) %>%
  arrange(desc(end.x)) %>%
  dplyr::mutate(rown = row_number())%>%
  dplyr::mutate(category = factor(category, levels = c("full-length", 
                                                       "5_degraded", 
                                                       "3_degraded", 
                                                       "degraded")))
  
pdf(here("figures/Fig5B_hsp20_read_classes.pdf"), height = 4,width = 4)
p <- ggplot(data = total_frame_sites_f2_plot %>% 
         mutate(condition = factor(condition, levels = c("WT", "dRNaseR"))),
       aes(color = condition)) +
  facet_grid(cols = vars(condition), scales = "free",
             rows = vars(category)) +
  geom_segment(aes(y = rown, yend = rown, x = start.x, xend = end.x),
               alpha = 1, size = 1) +
  scale_color_manual(values = c("#A9A9A9", "#D55E00")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_void()
print(p)
dev.off()
