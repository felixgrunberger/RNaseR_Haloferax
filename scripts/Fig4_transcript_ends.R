# Figure 4: Transcript end analysis

library(here)
source(here("scripts/helpers.R"))
dir.create(here("figures"), showWarnings = FALSE)

# Figure 4B and 4D — overlap of primary 5' and 3' ends
for (end in c("5", "3")) {
  peak_file <- if (end == "5") "tss_peaks_all_longread.tsv" else "tts_peaks_all_longread_cov5.tsv"
  class_col <- if (end == "5") "TSS_class" else "TTS_class"
  peaks <- vroom::vroom(here("data/Rdata", peak_file), show_col_types = FALSE) %>%
    dplyr::filter(biotype == "mRNA ORF", .data[[class_col]] == "primary") %>%
    dplyr::distinct(gene_id, set)
  stopifnot(all(c("twt", "tRNaseR") %in% peaks$set))
  sets <- list(WT = peaks$gene_id[peaks$set == "twt"],
               delta_rnr = peaks$gene_id[peaks$set == "tRNaseR"])
  p <- ggvenn::ggvenn(sets, fill_color = c("grey20", "#B7AEE3"),
                     stroke_size = 0.5, text_size = 5)
  ggplot2::ggsave(here("figures", paste0("Fig4", if (end == "5") "B" else "D", "_overlap.pdf")),
                  p, width = 4, height = 4)
}

# Figure 4C and 4E — positional distributions of terminal read counts.
# Source logic: 250131_end5_coverage.R / 250131_end3_coverage.R.
get_cov_bam <- function(plus_file, minus_file, barcode) {
  cov_one <- function(filename, strand) {
    x <- vroom::vroom(filename, col_names = FALSE, show_col_types = FALSE) %>%
      dplyr::rename(seqnames = X1, pos1 = X2, pos2 = X3, cov = X4) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(pos = list((pos1 + 1):pos2)) %>%
      tidyr::unnest(pos) %>%
      dplyr::ungroup()
    shifts <- if (strand == "forward") {
      list(cov_5down = dplyr::lag(x$cov, 4), cov_4down = dplyr::lag(x$cov, 3),
           cov_3down = dplyr::lag(x$cov, 2), cov_2down = dplyr::lag(x$cov, 1),
           cov_1up = dplyr::lead(x$cov, 1), cov_2up = dplyr::lead(x$cov, 2),
           cov_3up = dplyr::lead(x$cov, 3), cov_4up = dplyr::lead(x$cov, 4))
    } else {
      list(cov_5down = dplyr::lead(x$cov, 4), cov_4down = dplyr::lead(x$cov, 3),
           cov_3down = dplyr::lead(x$cov, 2), cov_2down = dplyr::lead(x$cov, 1),
           cov_1up = dplyr::lag(x$cov, 1), cov_2up = dplyr::lag(x$cov, 2),
           cov_3up = dplyr::lag(x$cov, 3), cov_4up = dplyr::lag(x$cov, 4))
    }
    dplyr::bind_cols(x, tibble::as_tibble(shifts)) %>%
      dplyr::mutate(strand = strand)
  }
  dplyr::bind_rows(cov_one(plus_file, "forward"), cov_one(minus_file, "reverse")) %>%
    dplyr::mutate(set = barcode)
}

barcode_conditions <- tibble::tibble(set = paste0("barcode", 20:23),
                                      cond = c("wt", "wt", "d", "d"))

for (end in c("5", "3")) {
  peak_file <- if (end == "5") "tss_peaks_all_longread.tsv" else "tts_peaks_all_longread_cov5.tsv"
  name_col <- if (end == "5") "TSS" else "TTS"
  class_col <- if (end == "5") "TSS_class" else "TTS_class"
  select_class <- if (end == "5") "primary" else "intragenic"
  peaks_all <- vroom::vroom(here("data/Rdata", peak_file), show_col_types = FALSE) %>%
    dplyr::filter(biotype == "mRNA ORF", .data[[class_col]] == select_class)
  common_gene_ids <- peaks_all %>%
    dplyr::distinct(gene_id, set) %>%
    dplyr::count(gene_id) %>%
    dplyr::filter(n == 2) %>%
    dplyr::pull(gene_id)
  peaks <- peaks_all %>%
    dplyr::filter(set == "twt", gene_id %in% common_gene_ids) %>%
    dplyr::distinct(.data[[name_col]], cov_peak, .keep_all = TRUE)

  bgs <- list.files(here("data/bamCoverage"),
                    pattern = paste0("end",end,".sorted\\.(pos|neg)\\.cpm\\.bedgraph$"),
                    full.names = TRUE, recursive = TRUE)
  full_cov <- purrr::map_dfr(barcode_conditions$set, function(bc) {
    pos <- bgs[stringr::str_detect(bgs, bc) & stringr::str_detect(bgs, "\\.pos\\.")]
    neg <- bgs[stringr::str_detect(bgs, bc) & stringr::str_detect(bgs, "\\.neg\\.")]
    if (length(pos) != 1 || length(neg) != 1)
      stop("Expected exactly one positive and one negative bedGraph for ", bc, "/end", end)
    get_cov_bam(pos, neg, bc)
  })

  peak_sites <- peaks %>%
    dplyr::select(gene_id, seqnames, strand, dplyr::all_of(name_col)) %>%
    dplyr::distinct()
  site_cov <- purrr::map_dfr(barcode_conditions$set, function(bc) {
    cov <- full_cov %>% dplyr::filter(set == bc)
    plus <- peak_sites %>% dplyr::filter(strand == "+") %>%
      dplyr::left_join(cov %>% dplyr::filter(strand == "forward"),
                       by = stats::setNames(c("seqnames", "pos"), c("seqnames", name_col)))
    minus <- peak_sites %>% dplyr::filter(strand == "-") %>%
      dplyr::left_join(cov %>% dplyr::filter(strand == "reverse"),
                       by = stats::setNames(c("seqnames", "pos"), c("seqnames", name_col)))
    dplyr::bind_rows(plus, minus) %>% dplyr::mutate(set = bc)
  })
  counts <- site_cov %>%
    dplyr::select(gene_id, set, cov, cov_5down, cov_4down, cov_3down, cov_2down,
                  cov_1up, cov_2up, cov_3up, cov_4up) %>%
    tidyr::pivot_longer(cols = c(cov, cov_5down, cov_4down, cov_3down,
                                 cov_2down, cov_1up, cov_2up, cov_3up, cov_4up),
                        names_to = "position_name", values_to = "counts") %>%
    dplyr::mutate(position = dplyr::recode(position_name,
      cov_5down = -4L, cov_4down = -3L, cov_3down = -2L, cov_2down = -1L,
      cov = 0L, cov_1up = 1L, cov_2up = 2L, cov_3up = 3L, cov_4up = 4L)) %>%
    dplyr::group_by(gene_id, set) %>%
    dplyr::mutate(rel_counts = counts / sum(counts, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(barcode_conditions, by = "set") %>%
    dplyr::group_by(gene_id, cond, position) %>%
    dplyr::summarise(mean_frac = mean(rel_counts, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(counts,
           ggplot2::aes(x = factor(position), y = mean_frac, fill = cond,
                        group = interaction(cond, position))) +
    ggplot2::stat_summary(alpha = 0.3, fun = median, geom = "bar",
                          position = ggplot2::position_dodge(), width = 0.8) +
    ggplot2::geom_boxplot(width = 0.8,
                          position = ggplot2::position_dodge(width = 0.8), outlier.colour = NA) +
    ggpubr::stat_compare_means(ggplot2::aes(group = cond), method = "wilcox.test", label = "p.signif") +
    ggplot2::scale_fill_manual(values = c(wt = "#A9A9A9", d = "#D55E00")) +
    ggpubr::theme_pubclean()
  ggplot2::ggsave(here("figures", paste0("Fig4", if (end == "5") "C" else "E", "_end_positions.pdf")),
                  p, width = 5, height = 6)
}
