library(tidyverse)

setwd("C:/Users/agb214/R files/ARG Project/data/contigs/")

F1A = read_tsv("F1A_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>%
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F1A")
F1B = read_tsv("F1B_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F1B")
F1C = read_tsv("F1C_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F1C")
F1D = read_tsv("F1D_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F1D")
F2A = read_tsv("F2A_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F2A")
F2B = read_tsv("F2B_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F2B")
F2C = read_tsv("F2C_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F2C")
F2D = read_tsv("F2D_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F2D")
F3A = read_tsv("F3A_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F3A")
F3B = read_tsv("F3B_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F3B")
F3C = read_tsv("F3C_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F3C")
F3D = read_tsv("F3D_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F3D")
F4A = read_tsv("F4A_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F4A")
F4B = read_tsv("F4B_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F4B")
F4C = read_tsv("F4C_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F4C")
F4D = read_tsv("F4D_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F4D")
F6A = read_tsv("F6A_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F6A")
F6B = read_tsv("F6B_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F6B")
F6C = read_tsv("F6C_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F6C")
F6D = read_tsv("F6D_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F6D")
F8A = read_tsv("F8A_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F8A")
F8B = read_tsv("F8B_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F8B")
F8C = read_tsv("F8C_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F8C")
F8D = read_tsv("F8D_taxa_TPM_contigs.tsv", na = c("","NA")) %>% 
  filter(TPM_abundance > 0) %>% 
  select(Contig, TPM_abundance, lineage) %>% mutate(sample = "F8D")

# define the levels of taxa in "Taxonomy"
taxa = c("superkingdom","phylum","class","order","family","genus","species")

contigs = rbind(F1A,F1B,F1C,F1D,
            F2A,F2B,F2C,F2D,
            F3A,F3B,F3C,F3D,
            F4A,F4B,F4C,F4D,
            F6A,F6B,F6C,F6D,
            F8A,F8B,F8C,F8D) %>%
  separate(lineage, taxa, ";") %>% # separate the 7 levels of taxa from lineage as new columns
  replace_na(setNames(as.list(rep("unclassified", length(taxa))), taxa)) %>% # remove any "NA" in taxa columns with "Unclassified"
  mutate(pond = substr(sample, 1,2),
         superkingdom = str_replace_all(superkingdom, "unclassified cellular organisms superkingdom","unclassified"),
         superkingdom = str_replace_all(superkingdom, "unclassified root superkingdom","unclassified"))

genus_count = contigs %>% select(Contig, TPM_abundance, superkingdom, genus, sample) %>%
  group_by(sample, superkingdom, genus) %>%
  summarise(count = mean(TPM_abundance), .groups = "drop") %>% # group by genus and 
  group_by(sample)

ctg_rel_abund_genus = genus_count %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count)
