library(tidyverse)

setwd("C:/Users/agb214/R files/ARG Project/data/viruses/")

F1A = read.csv("F1A_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F1A")
F1B = read.csv("F1B_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F1B")
F1C = read.csv("F1C_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F1C")
F1D = read.csv("F1D_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F1D")
F2A = read.csv("F2A_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F2A")
F2B = read.csv("F2B_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F2B")
F2C = read.csv("F2C_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F2C")
F2D = read.csv("F2D_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F2D")
F3A = read.csv("F3A_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F3A")
F3B = read.csv("F3B_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F3B")
F3C = read.csv("F3C_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F3C")
F3D = read.csv("F3D_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F3D")
F4A = read.csv("F4A_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F4A")
F4B = read.csv("F4B_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F4B")
F4C = read.csv("F4C_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F4C")
F4D = read.csv("F4D_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F4D")
F6A = read.csv("F6A_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F6A")
F6B = read.csv("F6B_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F6B")
F6C = read.csv("F6C_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F6C")
F6D = read.csv("F6D_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F6D")
F8A = read.csv("F8A_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F8A")
F8B = read.csv("F8B_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F8B")
F8C = read.csv("F8C_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F8C")
F8D = read.csv("F8D_RGI_BLAST.txt", sep="\t", na = c("","NA")) %>% 
  select("ORF_ID","Contig","Cut_Off","Best_Hit_ARO","ARO","Drug.Class",
         "Nudged") %>% mutate(sample = "F8D")

plasmids = rbind(F1A,F1B,F1C,F1D,
                 F2A,F2B,F2C,F2D,
                 F3A,F3B,F3C,F3D,
                 F4A,F4B,F4C,F4D,
                 F6A,F6B,F6C,F6D,
                 F8A,F8B,F8C,F8D)

viruses = plasmids[plasmids$Cut_Off %in% "Strict", ]

num_vir = length(unique(gsub("\\|.*","",plasmids$Contig)))

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