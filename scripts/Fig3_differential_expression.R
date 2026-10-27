# Figure 3: Differential gene expression

library(here)
source(here("scripts/helpers.R"))

# data ----
## gff ====
hvo_fasta <- readDNAStringSet(here("data/genome/hvo.fasta"))
names(hvo_fasta) <- str_split_fixed(names(hvo_fasta), " ",2)[,1]

### like TBA ####
hvo_gff <- read_in_gff_hvo(here("data/genome/hvo.gff"))

custom_gtf <- hvo_gff %>%
  dplyr::mutate(gene_id = locus_tag) %>%
  dplyr::rename("1" = seqid.x,
                "3" = type, 
                "4" = start, 
                "5" = end,
                "7" = strand) %>%
  dplyr::mutate("2" = "Genbank",
                "6" = "NA", 
                "8" = "NA") %>%
  dplyr::mutate("9" = paste0("transcript_id ",gene_id, "; ",
                             "gene_id ",gene_id, "; ",
                             "gene_name ",gene_id, ";")) %>%
  dplyr::select(-gene_id) %>%
  dplyr::select("1", "2", "3", "4", "5", "6", "7", "8", "9")

write.table(custom_gtf, file = here("data/genome/hvo_custom_all.gtf"),
            col.names = F, row.names = F, quote = F, sep = "\t")

gtf_file <- here("data/genome/hvo_custom_all.gtf")

# DESeq2 pipeline ----
dir <- paste0(here("data"), "/")


## mapped files ====
bam_list    <- list.files(path = here("data/mapped_fl_dedup"), 
                          recursive = T, full.names = T, pattern = ".sorted.bam$")
bam_list_f <- bam_list[!str_detect(bam_list, c("remap"))]
bam_list_f <- bam_list_f[!str_detect(bam_list_f, c("q20"))][1:12]

## calculate count matrix using featurecounts and label features ====
### calc ####
counts_all <- featureCounts(bam_list_f, verbose = F, annot.ext = gtf_file, strandSpecific = 1,
                            primaryOnly = T,isGTFAnnotationFile = T, nthreads = 8,
                            GTF.featureType = c("CDS","SRP_RNA", "RNase_P_RNA", "tRNA", "rRNA"), 
                            allowMultiOverlap = T, 
                            isLongRead = T)$counts

colnames(counts_all) <- str_remove_all(str_split_fixed(colnames(counts_all), "\\.",4)[,1], "_umi_dedup")

### normalise for export #####
counts_all_tab <- counts_all %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "barcode", values_to = "counts") %>%
  left_join(hvo_gff %>% dplyr::mutate(width = abs(start-end)), by = c("gene"="locus_tag")) %>%
  group_by(barcode) %>%
  mutate(rpk = (counts/width*1000),
         TPM_hand = rpk/(sum(rpk, na.rm = T)/1000000),
         rpkm = (counts/(sum(counts)/1000000))/width*1000) %>%
  dplyr::select(gene, barcode, type,counts,rpkm, TPM_hand)


vroom_write(counts_all_tab, here("data/Rdata/counts_all_tab.tsv"))
counts_all_tab <- vroom(here("data/Rdata/counts_all_tab.tsv"))

### merge 

### samples ####
ont_run_info <- read_xlsx(here("data/samples/ont_run_info.xlsx"))

# mRNA analysis ----
salmon_samples <- ont_run_info %>%
  dplyr::mutate(condition = as.factor(condition)) %>%
  dplyr::select(barcode, condition, replicate) %>%
  dplyr::filter(!str_detect(condition, "total"))
remove_rna_samples <- data.table(condition = c("WT",
                                               "dRNaseR"),
                                 replicate = c("1","3"),
                                 remove = T) %>%
  left_join(salmon_samples) %>%
  dplyr::select(barcode) %>%
  deframe()

salmon_samples_filter <- salmon_samples %>%
  dplyr::filter(!barcode %in% remove_rna_samples) 

## with replicate filter ====
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_all[,colnames(counts_all) %in% salmon_samples_filter$barcode], 
                                      colData = salmon_samples_filter,
                                      design = ~condition)

dds <- dds[rowSums(counts(dds)) > 1,]
## count table ====
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds,
                    coef = "condition_WT_vs_dRNaseR", 
                    type="ashr")

deseq_wt_rnaseR <- resLFC %>%
  as.data.frame() %>%
  rownames_to_column("locus_tag") %>%
  dplyr::mutate(log2FoldChange = -log2FoldChange) %>%
  left_join(hvo_gff, by = c("locus_tag" = "locus_tag"))

vroom_write(deseq_wt_rnaseR, paste0(dir, "Rdata/deseq_wt_rnaseR.tsv"))

deseq_wt_rnaseR <- vroom(paste0(dir, "Rdata/deseq_wt_rnaseR.tsv"))

deseq_wt_rnaseR_f <- deseq_wt_rnaseR %>%
  dplyr::mutate(reg_type = ifelse(padj < 0.05 & log2FoldChange >= 2, "up",
                                  ifelse(padj < 0.05 & log2FoldChange <= -2, "down","rest")))
hvo_annotation_arcog <- vroom(here("data/genome/hvo_annotation_arcog.tsv"))

dir.create(here("figures"), showWarnings = FALSE)
pdf(here("figures/Fig3_DESeq_MA.pdf"), height = 4,width = 6)
p <- ggplot(data = deseq_wt_rnaseR_f %>%
              dplyr::filter(type %in% "CDS") %>% 
         left_join(hvo_annotation_arcog %>% distinct(shortarcog, old_locus_tag,arcog_cat_big_name), by = c("old_locus_tag")),
       aes(x = log10(baseMean+1), y = log2FoldChange, color = reg_type, 
           text = paste0(locus_tag, "_", product))) +
  facet_grid(cols = vars(type)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2.5) +
  geom_point(data = deseq_wt_rnaseR_f %>% 
               left_join(hvo_annotation_arcog %>% distinct(shortarcog, old_locus_tag,arcog_cat_big_name), by = c("old_locus_tag")) %>%
               dplyr::filter(type %in% c("CDS"))%>%
               dplyr::filter(reg_type %in% c("down", "up")), 
             shape = 21, aes(fill = arcog_cat_big_name), color = "black", size = 4) +
  theme_minimal() +
  scale_color_manual(values = c( "#F8A280", "grey70","#1C5560", "grey80"), na.value = "grey80") +
  scale_fill_manual(values = c("#FFFFFF","#41412F", "#999999", "grey60"), na.value = "grey60") +
  theme(panel.grid.minor = element_blank()) +
  guides(color = "none") +
  theme_pubclean()
print(p)
dev.off()
