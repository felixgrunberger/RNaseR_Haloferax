# R packages and annotation functions

library(here)
library(tidyverse)
library(vroom)
library(data.table)
library(readxl)
library(ape)
library(Biostrings)
library(Rsubread)
library(DESeq2)
library(ashr)
library(GenomicAlignments)
library(Rsamtools)
library(ggpubr)

read_in_gff_hvo <-  function(input_gff){ 
  
  gff_p <- read.gff(input_gff) %>%
    dplyr::filter(type %in% c("CDS", "tRNA", "snoRNA", "rRNA", "RNase_P_RNA", "SRP_RNA")) %>%
    as_tibble() %>%
    dplyr::mutate(start_product = start, end_product = end) %>%
    dplyr::mutate(gene = str_split_fixed(str_split_fixed(attributes, ";Parent=",2)[,2], ";Dbxref",2)[,1],
                  product = str_split_fixed(str_split_fixed(attributes, "product=",2)[,2], ";",2)[,1]) %>%
    dplyr::select(seqid, gene, product, start, end, strand, type)
  
  gff_g <- read.gff(input_gff) %>%
    dplyr::filter(type %in% c("pseudogene","gene")) %>%
    as_tibble() %>%
    dplyr::mutate(gene = str_split_fixed(str_split_fixed(attributes, "ID=",2)[,2], ";Dbxref",2)[,1],
                  new_id = str_split_fixed(str_split_fixed(attributes, "Name==",2)[,2], ";gbkey",2)[,1], 
                  locus_tag = str_split_fixed(str_split_fixed(attributes, "locus_tag=",2)[,2], ";",2)[,1],
                  old_locus_tag = str_split_fixed(attributes, "old_locus_tag=",2)[,2]) %>%
    dplyr::select(seqid, gene, new_id, locus_tag, old_locus_tag)
  
  gff_final <- left_join(gff_g, gff_p, by = "gene") %>%
    dplyr::distinct(gene, .keep_all = T)
  
  
  return(gff_final)
}


read_in_gff_rrna <- function(input_file) {
  # Only the reference boundaries are used in the Figure 8 classification.
  ape::read.gff(input_file) %>%
    dplyr::filter(type == "rRNA") %>%
    dplyr::mutate(start_feature = start, end_feature = end) %>%
    dplyr::arrange(start_feature)
}
