library(tidyverse)

setwd("C:/Users/agb214/R files/ARG Project/data/plasmids/")

F1A = read.csv("F1A_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F1A")
F1B = read.csv("F1B_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F1B")
F1C = read.csv("F1C_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F1C")
F1D = read.csv("F1D_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F1D")
F2A = read.csv("F2A_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F2A")
F2B = read.csv("F2B_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F2B")
F2C = read.csv("F2C_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F2C")
F2D = read.csv("F2D_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F2D")
F3A = read.csv("F3A_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F3A")
F3B = read.csv("F3B_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F3B")
F3C = read.csv("F3C_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F3C")
F3D = read.csv("F3D_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F3D")
F4A = read.csv("F4A_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F4A")
F4B = read.csv("F4B_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F4B")
F4C = read.csv("F4C_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F4C")
F4D = read.csv("F4D_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F4D")
F6A = read.csv("F6A_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F6A")
F6B = read.csv("F6B_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F6B")
F6C = read.csv("F6C_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F6C")
F6D = read.csv("F6D_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F6D")
F8A = read.csv("F8A_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F8A")
F8B = read.csv("F8B_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F8B")
F8C = read.csv("F8C_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F8C")
F8D = read.csv("F8D_output_trm.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged","lineage","TPM_abundance") %>% mutate(sample = "F8D")

# define the levels of taxa in "Taxonomy"
taxa = c("superkingdom","phylum","class","order","family","genus","species")

plasmids = rbind(F1A,F1B,F1C,F1D,
                 F2A,F2B,F2C,F2D,
                 F3A,F3B,F3C,F3D,
                 F4A,F4B,F4C,F4D,
                 F6A,F6B,F6C,F6D,
                 F8A,F8B,F8C,F8D) %>% 
  separate(lineage, taxa, ";")

# remove any "NA" in taxa columns with "Unclassified"
plasmids[taxa][is.na(plasmids[taxa])] = "Unclassified"

# remove plasmids with viruses as taxnomic identity
plasmids = plasmids[! plasmids$superkingdom %in% "Viruses", ]


#remove the 1 plasmids with Nudged == True
plasmids_ARGs = plasmids[!plasmids$Nudged %in% "True",]

# remove plasmids with no abundance
plasmids_ARGs = plasmids_ARGs[plasmids_ARGs$TPM_abundance > 0, ]

# remove loose cut-off
plas_strict = plasmids_ARGs[! plasmids_ARGs$Cut_Off %in% "Loose",]

# if a gene conferes multiple resistances, split into new row copies with each resistance
plas_strict = plas_strict %>% mutate(Drug.Class = strsplit(as.character(Drug.Class), "; ")) %>% 
  unnest(Drug.Class)

# remove the word "antibiotic" and any spaces from the Drug.Class column
plas_strict$Drug.Class = str_trim(gsub("antibiotic", "", plas_strict$Drug.Class), side = c("both"))

# remove NAs
plas_strict = plas_strict[!is.na(plas_strict$species), ]
