# Figure 8: rRNA precursor classification

library(here)
source(here("scripts/helpers.R"))

# functions ----
read_in_bam_large <- function(file){
  
  init <- readGAlignments(file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
  init_t <- GenomicAlignments::as.data.frame(init) %>%
    dplyr::mutate(minion_read_name = names(init),
                  mapped_gene = seqnames) 
  
  left  <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  init_t$soft_l <- as_tibble(cigarOpTable(left))$S
  init_t$hard_l <- as_tibble(cigarOpTable(left))$H
  init_t$soft_r <- as_tibble(cigarOpTable(right))$S
  init_t$hard_r <- as_tibble(cigarOpTable(right))$H
  init_t$matches <- as_tibble(cigarOpTable(init_t$cigar))$M
  init_t$insertions <- as_tibble(cigarOpTable(init_t$cigar))$I
  init_t$deletions <- as_tibble(cigarOpTable(init_t$cigar))$D
  init_t$skip <- as_tibble(cigarOpTable(init_t$cigar))$N
  init_t$cigar_S <- as_tibble(cigarOpTable(init_t$cigar))$S
  
  # calculate number of aligned reads based on CIGAR operations (M,I)
  init_t$aligned_reads <- unlist(lapply(explodeCigarOpLengths(init_t$cigar, ops = c("M", "I")), function(x) sum(x)))
  
  # calc read identity..
  init_t_final <- init_t %>%
    dplyr::mutate(identity = (1 - NM/aligned_reads)*100) %>%
    dplyr::group_by(minion_read_name) %>%
    dplyr::filter(identity == max(identity),
                  aligned_reads == max(aligned_reads)) %>%
    dplyr::distinct(minion_read_name, .keep_all = T) %>%
    dplyr::mutate(gene = str_split_fixed(mapped_gene,"-",2)[,2])
  
  return(init_t_final)
}


group_maps16 <- function(input,bc){
  
  border <- 500
  m1 <-  border - 6 
  m2 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2
  m3 <- (1599741)-(hvo_gff_rrna$end_feature[1]-border)+2+(hvo_gff_rrna$start_feature[1]-1598084)
  c1 <- 20
  
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= 90, strand == "+") %>%
    dplyr::mutate(group = case_when((start >= m2 & start < m3 & end >= m3 & soft_l >= c1) ~ "(01) 5_extended_16S",
                                    (start >= m2 & start < m3 & end >= m3 & soft_l < c1) ~ "(02) post-16S-bhb",
                                    (start <= m1 & end >= m3) ~ "(03) circ_16S",
                                    (start <= m1 & end >= m2 & end < m3 & soft_r < c1) ~ "(04) open_circ_16S_5",
                                    (end <= m1 & soft_r <= c1 & soft_l < 2000) ~ "(06) mature_16S",
                                    (start >= m3 & soft_l <= c1 & soft_r < 2000) ~ "(06) mature_16S")) %>%
    dplyr::mutate(barcode = bc)
}

group_maps23 <- function(input,bc){
  
  border <- 500
  m1 <- border + 33 
  m2 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2
  m3 <- (hvo_gff_rrna$end_feature[2]+72)-(hvo_gff_rrna$end_feature[2]-border)+2+(hvo_gff_rrna$start_feature[2]-11-(hvo_gff_rrna$start_feature[2]-50))
  c1 <- 20
  
  mat_l <- (hvo_gff_rrna$end_feature[2]+33)-(hvo_gff_rrna$start_feature[2]-11)
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= 90, strand == "+") %>%
    dplyr::mutate(group = case_when((start < m1 & end > m3) ~ "(10) circ_23S",
                                    (start < (m1) & end > m1 & end <= m2 & soft_r <= c1) ~ "(11) open_circ_23S_cutBHB3",
                                    (start < m1  & end <= m1 & soft_r <= c1 & width < mat_l) ~ "(12) mature_23S",
                                    (start >= m3  & end > m3 & soft_l <= c1 & width < mat_l) ~ "(12) mature_23S")) %>%
    dplyr::mutate(barcode = bc)
}

group_maps_all <- function(input,bc){
  
  bhb1 <- 1598084
  bhb2 <- 1599741
  bhb3 <- hvo_gff_rrna$start_feature[2] - 50 - 11
  bhb4 <- hvo_gff_rrna$end_feature[2] + 72
  mat1 <- hvo_gff_rrna$start_feature[1]
  mat2 <- hvo_gff_rrna$end_feature[1]-8
  mat3 <- hvo_gff_rrna$start_feature[2] - 11
  mat4 <- hvo_gff_rrna$end_feature[2] + 32
  mat5 <- hvo_gff_rrna$start_feature[3]
  mat6 <- hvo_gff_rrna$end_feature[3]
  trna_end <- 1599833
  c1 <- 20
  
  read_in_bam_large(input) %>%
    ungroup() %>%
    dplyr::filter(identity >= 90, strand == "+") %>%
    dplyr::mutate(group = case_when((start < bhb1 & start >= (mat1-444) & end > mat1 & soft_l <= c1 & soft_r <= c1) ~ "(01) 5_extended_16S",
                                    (start >= bhb1 & start < mat1 & end > mat1 & end <= bhb2 & soft_l <= c1 & soft_r <= c1) ~ "(02) post-16S-bhb",
                                    (start >= mat1 & end > mat2 & end <= bhb2 & soft_r < c1 & soft_l < c1) ~ "(05) open_circ_16S_cutBHB3",
                                    (start >= mat1 & start < mat2 & end > mat1 & end <= mat2 & soft_l <= c1 & soft_r <= c1) ~ "(06) mature_16S",
                                    (start >= trna_end & start < bhb3 & end > mat3 & soft_l <= c1 & soft_r <= c1) ~ "(07) 5_extended_23S",
                                    (start >= trna_end & end > bhb4 & end <= mat5 & soft_l <= c1 & soft_r <= c1) ~ "(08) 3_extended_23S",
                                    (start >= bhb3 & start < mat3 & end > mat4 & end <= bhb4 & soft_l < c1 & soft_r < c1) ~ "(09) post-23S-bhb",
                                    (start >= mat3 & start < mat4 & end > mat4 & end <= bhb4 & soft_l < c1 & soft_r < c1) ~ "(11) open_circ_23S_cutBHB3",
                                    (start >= mat3 & end <= mat4 & soft_l <= c1 & soft_r <= c1) ~ "(12) mature_23S")) %>%
    dplyr::mutate(barcode = bc)
}


# data ----

## genome ====
hvo_gff_rrna <- read_in_gff_rrna(here("data/genome/hvo.gff"))

## circ category table ====
dir <- paste0(here("data"), "/")

### 16S ####
list16circ <- list.files(paste0(dir,"mapped_fl_only_umi_removed_circ16S"), recursive = T, pattern = ".sorted.bam$",full.names = T)
list16barc <- str_extract(list16circ, "barcode\\d+")
circ_map16_cat <- pmap_dfr(list(list16circ,
                                list16barc),
                           group_maps16)

### 23S ####
list23circ <- list.files(paste0(dir,"mapped_fl_only_umi_removed_circ23S"), recursive = T, pattern = ".sorted.bam$",full.names = T)
list23barc <- str_extract(list23circ, "barcode\\d+")
circ_map23_cat <- pmap_dfr(list(list23circ,
                                list23barc),
                           group_maps23)

## normal mapping table ====
listallcirc <- list.files(paste0(dir,"mapped_fl_only_umi_removed"), recursive = T, pattern = ".sorted.bam$",full.names = T)
listallcirc_f <- listallcirc[!str_detect(listallcirc, c("remap"))]
listallcirc_f <- listallcirc_f[!str_detect(listallcirc_f, c("q20"))][1:12]
listallbarc <- str_extract(listallcirc_f, "barcode\\d+")
all_map_cat <- pmap_dfr(list(listallcirc_f,
                             listallbarc),
                        group_maps_all)

## combine tables ====
comb_map <- all_map_cat %>%
  full_join(circ_map16_cat %>%
              dplyr::rename(group_16c = group) %>%
              dplyr::select(barcode,minion_read_name, group_16c),
            by = "minion_read_name") %>%
  full_join(circ_map23_cat %>%
              dplyr::rename(group_23c = group) %>%
              dplyr::select(barcode,minion_read_name, group_23c),
            by = "minion_read_name") %>%
  dplyr::mutate(group_all = ifelse(!is.na(group_16c), group_16c, group_23c),
                group_all = ifelse(is.na(group_all),group, group_all),
                clip_g = case_when(((soft_l > soft_r) & (soft_l > 20)) ~ "left-clipped",
                                   ((soft_l < soft_r) & (soft_l > 20)) ~ "right-clipped",
                                   ((soft_r > soft_l) & (soft_r > 20)) ~ "right-clipped",
                                   ((soft_r < soft_l) & (soft_r > 20)) ~ "left-clipped",
                                   (soft_l <= 20) ~ "none", 
                                   (soft_r <= 20) ~ "none"),
                group_all = ifelse((group == "(05) open_circ_16S_cutBHB3" & group_all == "(06) mature_16S"), "(05) open_circ_16S_cutBHB3",
                                   ifelse((group_all == "(01) 5_extended_16S" & (soft_r > 20 | soft_l > 20)), NA,
                                          ifelse((group_all == "(02) post-16S-bhb" & (soft_r > 20 | soft_l > 20)), NA,
                                                 ifelse((group_all == "(06) mature_16S" & (soft_r > 20 | soft_l > 20)),NA,
                                                        ifelse((group_all == "(11) open_circ_23S_cutBHB3" & (clip_g != "none")),NA,
                                                               ifelse((group_all == "(12) mature_23S" & (clip_g != "none")),NA,group_all)))))))

vroom_write(comb_map, paste0(dir, "/Rdata/comb_map.tsv"))

# RNA precursor count matrix, DESeq2 (Figure 8B)
ont_run_info <- readxl::read_xlsx(here("data/samples/ont_run_info.xlsx"))
table_pre <- comb_map %>%
  dplyr::mutate(rRNA_group = ifelse(str_detect(group_all, "23S"), "23S",
                                    ifelse(str_detect(group_all, "16S"), "16S","else"))) %>%
  dplyr::mutate(bc = ifelse(!is.na(barcode.x),barcode.x,
                            ifelse(is.na(barcode.x) & !is.na(barcode.y), barcode.y, barcode))) %>%
  dplyr::filter(!is.na(bc),
                bc != "NA") %>%
  left_join(ont_run_info, by = c("bc" = "barcode")) %>%
  group_by(bc, group_all) %>%
  summarise(counts = n()) %>%
  dplyr::filter(!is.na(bc),
                bc != "NA") %>%
  distinct(bc, group_all, counts) %>%
  dplyr::filter(!is.na(bc),
                bc != "NA") %>%
  group_by(bc) %>%
  dplyr::mutate(frac = counts/sum(counts)*100) %>%
  drop_na(group_all) %>%
  dplyr::filter(bc %in% paste0("barcode", 13:18)) %>%
  #dplyr::mutate(condition = factor(condition, levels = c("WT_total", "dRNaseR_total"))) %>%
  dplyr::select(bc, counts, group_all) %>%
  pivot_wider(names_from = bc, values_from = counts)


vroom_write(table_pre, here("data/Rdata/table_pre.tsv"))
table_pre <- table_pre %>%
  drop_na()
table_pre_matrix <- as.matrix(table_pre[,-1])
rownames(table_pre_matrix) <- table_pre$group_all

table_info <- ont_run_info %>%
  dplyr::filter(barcode %in% colnames(table_pre_matrix)) %>%
  dplyr::select(barcode, condition) %>%
  tibble::column_to_rownames("barcode")
table_info <- table_info[colnames(table_pre_matrix), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = table_pre_matrix,
                              colData = table_info,
                              design= ~ condition)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef="condition_WT_total_vs_dRNaseR_total", type="ashr")
res_tab <- as.data.frame(res) %>%
  rownames_to_column("precursor") %>%
  dplyr::mutate(log2FoldChange = -log2FoldChange) %>%
  dplyr::mutate(pre = str_remove_all(str_remove_all(str_split_fixed(precursor, " ", 2)[,1], "\\("), "\\)")) %>%
  dplyr::mutate(reg_type = ifelse(padj < 0.05 & log2FoldChange > 0, "up",
                                  ifelse(padj < 0.05 & log2FoldChange < 0, "down","rest")))

dir.create(here("figures"), showWarnings = FALSE)
pdf(here("figures/Fig8B_rRNA_precursor_DESeq.pdf"), height = 5,width = 3.5)
p <- ggplot(data = res_tab,
       aes(x = log10(baseMean), y = log2FoldChange, fill = reg_type)) +
  geom_text(aes(label = pre), vjust = -1.5) +
  geom_hline(yintercept = 0, color = "black") +
  theme_pubclean() +
  geom_point(shape = 21, size = 5) +
  scale_fill_manual(values = c("grey65","#BFB8F1"), na.value = "grey75") +
  scale_y_continuous(limits = c(-1.5,1.5))
print(p)
dev.off()

# Transcript class percentages (Figure 8C)
counts_all_tab <- vroom::vroom(here("data/Rdata/counts_all_tab.tsv"), show_col_types = FALSE)
table_pre <- vroom(here("data/Rdata/table_pre.tsv"))
table_pre_g <- table_pre %>%
  pivot_longer(-1,names_to = "barcode", values_to = "counts") %>%
  dplyr::mutate(type = "rRNA_precursor",
                rpkm = NA,
                TPM_hand = NA) %>%
  dplyr::rename(gene = 1,
                barcode = 2) %>%
  dplyr::mutate(type = ifelse(gene %in% "(04) open_circ_16S_5", "rRNA_precursor_want", type)) %>%
  left_join(ont_run_info, by = "barcode")

merge_tab <- counts_all_tab %>%
  dplyr::filter(!gene %in% c("HVO_RS18920","HVO_RS13015")) %>%
  dplyr::mutate(barcode = str_remove_all(barcode, "_umi_dedup.sorted.bam")) %>%
  left_join(ont_run_info, by = "barcode")

table_pre_g2 <- table_pre_g %>%
  bind_rows(merge_tab) %>%
  group_by(condition) %>%
  dplyr::mutate(total = sum(counts)) %>%
  group_by(condition, type) %>%
  reframe(rel = sum(counts)/total*100) %>%
  distinct(condition, type, rel) %>%
  dplyr::mutate(condition = factor(condition, levels = c("WT", "dRNaseR",
                                                         "WT_total", "dRNaseR_total"))) %>%
  dplyr::mutate(type = factor(type, 
                              levels = rev(c("CDS", "RNase_P_RNA", "SRP_RNA",
                                         "rRNA", "rRNA_precursor", "rRNA_precursor_want"))))

pdf(here("figures/Fig8C_RNA_composition.pdf"), height = 6,width = 6)
p <- ggplot(data = table_pre_g2 %>%
         dplyr::filter(str_detect(condition, "total")),
       aes(x = condition, y = rel, fill = type)) +
  geom_bar(stat = "identity", color = "black") +
  theme_pubclean()
print(p)
dev.off()
