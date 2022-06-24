library(tidyverse)

########## ARGs ##########

setwd("C:/Users/agb214/R files/ARG Project/")

SAMPLE = read.csv("data/genecalls/SAMPLE_taxa_TPM_genecalls.tsv", sep = "\t", na.strings=c("","NA")) %>% 
  select("ORF_ID","Cut_Off","Best_Hit_ARO","ARO","Drug.Class","Nudged","lineage","TPM_abundance","TaxID") %>% 
  mutate(sample = "SAMPLE")

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

###
### preprocess read taxnonomy from kaiju and kraken output ###
###

taxa = c("superkingdom","phylum","class","order","family","genus","species")

SAMPLE = read_tsv("SAMPLE_combined_rmft.tsv", na = c("","NA"), col_names = c("count","lineage")) %>%
  count(lineage) %>%
  mutate(sample = "SAMPLE",
         lineage = str_replace_all(lineage, "unclassified cellular organisms","unclassified"),
         lineage = str_replace_all(lineage, "unclassified root","unclassified"),
         farm = substr(sample, 1,2),
         rel_abund = n / sum(n)) %>%
  separate(lineage, taxa, ";") %>% # separate the 7 levels of taxa from lineage as new columns
  replace_na(setNames(as.list(rep("unclassified", length(taxa))), taxa)) # remove any "NA" in taxa columns with "Unclassified"

df1 = rbind(F1A,F1B,F1C,F1D,
            F2A,F2B,F2C,F2D,
            F3A,F3B,F3C,F3D,
            F4A,F4B,F4C,F4D,
            F6A,F6B,F6C,F6D,
            F8A,F8B,F8C,F8D)
            
reads_by_genus = reads %>%
  group_by(sample, genus) %>%
  summarise(rel_abund = count/sum(count), .groups = "drop")


###
### preprocess plasmids from diamond identity 
### 

SAMPLE = read.csv("SAMPLE_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Best_Hit_ARO","Drug.Class","start","end","func","direction") %>% 
  mutate(sample = "SAMPLE")

plasmids = rbind(F1A,F1B,F1C,F1D,
                 F2A,F2B,F2C,F2D,
                 F3A,F3B,F3C,F3D,
                 F4A,F4B,F4C,F4D,
                 F6A,F6B,F6C,F6D,
                 F8A,F8B,F8C,F8D)
                 
# if a gene conferes multiple resistances, split into new row copies with each resistance
plas_strict = plasmids %>% 
  mutate(Drug.Class = strsplit(as.character(Drug.Class), "; ")) %>% 
  unnest(Drug.Class)

# remove the word "antibiotic" and any spaces from the Drug.Class column
plas_strict$Drug.Class = str_trim(gsub("antibiotic", "", plas_strict$Drug.Class), side = c("both"))
