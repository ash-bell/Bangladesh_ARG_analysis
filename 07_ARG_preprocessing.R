library(tidyverse)

########## ARGs ##########

setwd("C:/Users/agb214/R files/ARG Project/")

F1A = read.csv("data/genecalls/F1A_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F1A")
F1B = read.csv("data/genecalls/F1B_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F1B")
F1C = read.csv("data/genecalls/F1C_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F1C")
F1D = read.csv("data/genecalls/F1D_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F1D")
F2A = read.csv("data/genecalls/F2A_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F2A")
F2B = read.csv("data/genecalls/F2B_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F2B")
F2C = read.csv("data/genecalls/F2C_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F2C")
F2D = read.csv("data/genecalls/F2D_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F2D")
F3A = read.csv("data/genecalls/F3A_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F3A")
F3B = read.csv("data/genecalls/F3B_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F3B")
F3C = read.csv("data/genecalls/F3C_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F3C")
F3D = read.csv("data/genecalls/F3D_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F3D")
F4A = read.csv("data/genecalls/F4A_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F4A")
F4B = read.csv("data/genecalls/F4B_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F4B")
F4C = read.csv("data/genecalls/F4C_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F4C")
F4D = read.csv("data/genecalls/F4D_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F4D")
F6A = read.csv("data/genecalls/F6A_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F6A")
F6B = read.csv("data/genecalls/F6B_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F6B")
F6C = read.csv("data/genecalls/F6C_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F6C")
F6D = read.csv("data/genecalls/F6D_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F6D")
F8A = read.csv("data/genecalls/F8A_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F8A")
F8B = read.csv("data/genecalls/F8B_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F8B")
F8C = read.csv("data/genecalls/F8C_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F8C")
F8D = read.csv("data/genecalls/F8D_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "F8D")

# define the levels of taxa in "Taxonomy"
taxa = c("superkingdom","phylum","class","order","family","genus","species")

# combine df with taxa groupings as new columns
df1 = rbind(F1A,F1B,F1C,F1D,
            F2A,F2B,F2C,F2D,
            F3A,F3B,F3C,F3D,
            F4A,F4B,F4C,F4D,
            F6A,F6B,F6C,F6D,
            F8A,F8B,F8C,F8D) %>% 
  separate(lineage, taxa, ";") %>%
  .[!is.na(.$Drug.Class), ] # rows with all NAs get introduced. I don't know why or how but this removes them

# remove any "NA" in taxa columns with "Unclassified"
df1[taxa][is.na(df1[taxa])] = "Unclassified"

# remove no abundance
ARGs = ARGs[ARGs$TPM_abundance > 0, ]

### split multiple antibiotics ### 

# if a gene confers multiple resistances, split into new row copies with each resistance
ARGs.expand = ARGs %>% mutate(Drug.Class = strsplit(as.character(Drug.Class), "; ")) %>% unnest(Drug.Class)

# remove the word "antibiotic" and any spaces from the Drug.Class column
ARGs.expand$Drug.Class = str_trim(gsub("antibiotic", "", ARGs.expand$Drug.Class), side = c("both"))

### if split by multiple antibiotics also divide TPM by drugclass count (multiple resistances)
## therefore ARG TPM isn't doubled if split
#use if not plotting by drug class

ARGs_ARO_safe = ARGs %>% 
  mutate(count = str_count(Drug.Class, pattern = "; ")+1) %>%
  mutate(Drug.Class = strsplit(as.character(Drug.Class), "; ")) %>% 
  unnest(Drug.Class) %>%
  mutate(TPM_abundance = TPM_abundance/count)

ARGs_ARO_safe$Drug.Class = str_trim(gsub("antibiotic", "", ARGs_ARO_safe$Drug.Class), side = c("both"))
