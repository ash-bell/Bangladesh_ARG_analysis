### setup and preprocessing ###

library(tidyverse)
library(factoextra)
library(FactoMineR)
library(missMDA)
library(vegan)
library(cowplot)
library(metacoder)
library(ggrepel)
library(car)
library(rstatix)
library(lme4)
library(jtools)
library(pheatmap)
library(gggenes)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(sp)
library(maps)

### Setup ###
set.seed(1234)
setwd("C:/Users/agb214/R files/ARG Project/")
source("scripts/07_ARG_read_plasmid_preprocessing.R")

### Define my sample labels ###
sample_levels = c("F1_Mon","F1_PosMon","F1_Win","F1_PreMon",
                  "F2_Mon","F2_PosMon","F2_Win","F2_PreMon",
                  "F3_Mon","F3_PosMon","F3_Win","F3_PreMon",
                  "F4_Mon","F4_PosMon","F4_Win","F4_PreMon",
                  "F5_Mon","F5_PosMon","F5_Win","F5_PreMon",
                  "F6_Mon","F6_PosMon","F6_Win","F6_PreMon")

### Read Metadata table ###
metadata = read_csv("/metadata/BD_fish_ponds_pooled_samples_for_metagenome_seq.csv") %>%
  select("SampleID","Farmer","Ponds","Main_crop","Village","Upazila","Timepoints","Dates","Season","Samples_in_pool") %>%
  mutate(Season = gsub("_", " ", Season),
         Upazila = gsub("_", " ", Upazila)) %>% 
  mutate(SampleID = gsub("A$", "_Mon", SampleID)) %>%
  mutate(SampleID = gsub("B$", "_PosMon", SampleID)) %>%
  mutate(SampleID = gsub("C$", "_Win", SampleID)) %>%
  mutate(SampleID = gsub("D$", "_PreMon", SampleID)) %>%
  mutate(SampleID = gsub("^F8", "F5", SampleID)) %>%
  mutate(Farmer = gsub("^F8$", "F5", Farmer)) %>% # rename Farm F8 to F5 so farm order increases numerically
  arrange(SampleID)

metadata$SampleID = factor(metadata$SampleID, levels = sample_levels)

### Ensure the samples labels from all input data is defined how I want them ###
reads = reads %>%
  mutate(sample = gsub("A$", "_Mon", sample)) %>%
  mutate(sample = gsub("B$", "_PosMon", sample)) %>%
  mutate(sample = gsub("C$", "_Win", sample)) %>%
  mutate(sample = gsub("D$", "_PreMon", sample)) %>%
  mutate(sample = gsub("^F8", "F5", sample)) %>%
  mutate(farm = gsub("^F8$", "F5", farm)) %>%
  arrange(sample)

reads$sample = factor(reads$sample, levels = sample_levels)

ARGs.expand = ARGs.expand %>%
  mutate(sample = gsub("A$", "_Mon", sample)) %>%
  mutate(sample = gsub("B$", "_PosMon", sample)) %>%
  mutate(sample = gsub("C$", "_Win", sample)) %>%
  mutate(sample = gsub("D$", "_PreMon", sample)) %>%
  mutate(sample = gsub("^F8", "F5", sample)) %>%
  arrange(sample)

ARGs.expand$sample = factor(ARGs.expand$sample, levels = sample_levels)

ARGs_ARO_safe = ARGs_ARO_safe %>%
  mutate(sample = gsub("A$", "_Mon", sample)) %>%
  mutate(sample = gsub("B$", "_PosMon", sample)) %>%
  mutate(sample = gsub("C$", "_Win", sample)) %>%
  mutate(sample = gsub("D$", "_PreMon", sample)) %>%
  mutate(sample = gsub("^F8", "F5", sample)) %>%
  arrange(sample)

ARGs_ARO_safe$sample = factor(ARGs_ARO_safe$sample, levels = sample_levels)

reads_by_genus = reads_by_genus %>%
  mutate(sample = gsub("A$", "_Mon", sample)) %>%
  mutate(sample = gsub("B$", "_PosMon", sample)) %>%
  mutate(sample = gsub("C$", "_Win", sample)) %>%
  mutate(sample = gsub("D$", "_PreMon", sample)) %>%
  mutate(sample = gsub("^F8", "F5", sample)) %>%
  arrange(sample)

reads_by_genus$sample = factor(reads_by_genus$sample, levels = sample_levels)

KEGGs = KEGGs %>%
  mutate(SampleID = gsub("A$", "_Mon", SampleID)) %>%
  mutate(SampleID = gsub("B$", "_PosMon", SampleID)) %>%
  mutate(SampleID = gsub("C$", "_Win", SampleID)) %>%
  mutate(SampleID = gsub("D$", "_PreMon", SampleID)) %>%
  mutate(SampleID = gsub("^F8", "F5", SampleID)) %>%
  arrange(SampleID)

KEGGs$SampleID = factor(KEGGs$SampleID, levels = sample_levels)


### Define the colours of my villages and upazilas ### 
clrs = c(rep('#e6194B', 4),
         rep('#3cb44b', 4),
         rep('#ffe119', 4),
         rep('#4363d8', 4),
         rep('#911eb4', 4),
         rep('#f58231', 4))
names(clrs) = unique(reads$sample)
farms = unique(clrs)
names(farms) = c(1,2,3,4,5,6)
clrs = append(clrs,farms)
clrs$`Jamalpur Sadar` = "#e377c2"
clrs$Muktagacha = "#17becf"
clrs$Tarakanda = "#bcbd22"
colScale <- scale_colour_manual(name = "pond",values = clrs)
fillScale <- scale_fill_manual(name = "pond",values = clrs)

### Figure 1 
############ World Map #####################
world <- ne_countries(scale = "medium", returnclass = "sf")

# co-ordinates are rounded up to preseve farmer anonomity
metadata = 
  metadata %>%
  mutate(latitude = ifelse(Village == "Kandulia", 24.8,
                           ifelse(Village == "Malotipur", 24.7,
                                  ifelse(Village == "Bamunpara", 24.8,
                                         ifelse(Village == "Gopalpur", 24.5,
                                                ifelse(Village == "Fulbaria", 24.6, NA)))))) %>%
  mutate(longitude = ifelse(Village == "Kandulia", 90.7,
                          ifelse(Village == "Malotipur", 90.2,
                                 ifelse(Village == "Bamunpara", 90.0,
                                        ifelse(Village == "Gopalpur", 89.9,
                                               ifelse(Village == "Fulbaria", 90.2, NA)))))) %>%
  mutate(Village = ifelse(Village == "Kandulia", "Kandulia\n(F1, F2)",
                           ifelse(Village == "Malotipur", "Malatipur\n(F3)",
                                  ifelse(Village == "Bamunpara", "Bamunpara\n(F4)",
                                         ifelse(Village == "Gopalpur", "Gopalpur\n(F6)",
                                                ifelse(Village == "Fulbaria", "Fulbaria\n(F5)", NA))))))

village = distinct(metadata, Village, .keep_all = T)

metadata = st_as_sf(metadata, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

# Upazila Shapefile from here https://geodash.gov.bd/layers/geonode:administrative_boundary_of_bangladesh_upazila_level
Upazila2 = read_sf("~/../Downloads/administrative_boundary_of_bangladesh_upazila_level/administrative_boundary_of_bangladesh_upazila_level.shp")

p1 = 
ggplot(data = world) +
  geom_sf() +
  geom_sf(data = Upazila2, aes(fill = ifelse(NAME_4=="Muktagachha", 'Muktagacha',
                                             ifelse(NAME_4=="Jamalpur S.", 'Jamalpur Sadar',
                                                    ifelse(NAME_4=="Phulpur", 'Tarakanda',
                                                           ifelse(NAME_2=="Mymensingh", 'Mymensingh', 
                                                                  ifelse(NAME_2=="Jamalpur", 'Mymensingh', NA))))))) +
  scale_fill_manual(values=c('#e377c2', '#17becf', '#dddddd', '#bcbd22'), name="Regions") +
  coord_sf(xlim = c(88, 92.72), ylim = c(21, 26.5))

p2 = 
ggplot(data = world) +
  geom_sf() +
  geom_sf(data = Upazila2, aes(fill = ifelse(NAME_4=="Muktagachha", 'Muktagacha',
                                             ifelse(NAME_4=="Jamalpur S.", 'Jamalpur Sadar',
                                                    ifelse(NAME_4=="Phulpur", 'Tarakanda',
                                                           ifelse(NAME_2=="Mymensingh", 'Mymensingh', 
                                                                  ifelse(NAME_2=="Jamalpur", 'Mymensingh', NA))))))) +
  scale_fill_manual(values=c('#e377c2', '#17becf', '#dddddd', '#bcbd22'), name="Regions") +
  geom_sf(data = metadata, aes(size = 1)) +
  geom_label_repel(data = village, label.size = 2, aes(x=longitude, y=latitude, label=Village, colour=Village), force = 75) +
  scale_colour_manual(values=c('#4363d8', '#911eb4', '#f58231', '#e6194B', '#3cb44b', '#ffe119'), name="Villages") +
  coord_sf(xlim = c(89.7, 91.2), ylim = c(24.3,25.2)) +
  theme(legend.position = "bottom", legend.box = "vertical")

plot = plot_grid(p1, p2, labels = c("A","B"),
          rel_widths = c(1,3),
          rel_heights = c(3,1))


### Fish diagram (Figure 1B constructed in affinty Designer)

### Figure 2, Supplimentary Figures 1 and 2. ###
### Relative abundance barplots by ponds using reads ###

clrs2 = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
         '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', 
         '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
         '#000075', '#808080', '#ffffff', '#000000')

bac  = reads %>%
  select(superkingdom, phylum, sample, rel_abund) %>%
  filter(superkingdom == "Bacteria") %>%
  select(-superkingdom) %>%
  group_by(sample, phylum) %>%
  summarise(sum_rel_abund = sum(rel_abund), .groups = "drop_last") %>%
  mutate(perc_sum_rel_abund = sum_rel_abund/sum(sum_rel_abund)) %>%
  slice_max(sum_rel_abund, n=10) %>%
  group_modify(~ add_row(.x, phylum = "Other phyla", 
                         perc_sum_rel_abund = (1-sum(.$perc_sum_rel_abund)), 
                         .before = 0)) %>%
  ungroup() %>%
  select(-sum_rel_abund)

bac.plot = ggplot(bac, aes(x=as.factor(sample), y=perc_sum_rel_abund*100, fill=phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage relative abundance (reads)",
       x = "Sample",
       fill="Phylum") +
  scale_fill_manual(values = clrs2) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust = 0),
        legend.position = "bottom")

vir = reads %>%
  select(superkingdom, phylum, sample, rel_abund) %>%
  filter(superkingdom == "Viruses") %>%
  select(-superkingdom) %>%
  group_by(sample, phylum) %>%
  summarise(sum_rel_abund = sum(rel_abund), .groups = "drop_last") %>%
  mutate(perc_sum_rel_abund = sum_rel_abund/sum(sum_rel_abund)) %>%
  slice_max(sum_rel_abund, n=10) %>%
  group_modify(~ add_row(.x, phylum = "Other phyla", 
                         perc_sum_rel_abund = (1-sum(.$perc_sum_rel_abund)), 
                         .before = 0)) %>%
  ungroup() %>%
  select(-sum_rel_abund)

vir.plot = ggplot(vir, aes(x=sample, y=perc_sum_rel_abund*100, fill=phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage relative abundance (reads)",
       x = "Sample",
       fill="Phylum") +
  scale_fill_manual(values = clrs2) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust = 0),
        legend.position = "bottom")

arc = reads %>%
  select(superkingdom, phylum, sample, rel_abund) %>%
  filter(superkingdom == "Archaea") %>%
  select(-superkingdom) %>%
  group_by(sample, phylum) %>%
  summarise(sum_rel_abund = sum(rel_abund), .groups = "drop_last") %>%
  mutate(perc_sum_rel_abund = sum_rel_abund/sum(sum_rel_abund)) %>%
  slice_max(sum_rel_abund, n=10) %>%
  group_modify(~ add_row(.x, phylum = "Other phyla", 
                         perc_sum_rel_abund = (1-sum(.$perc_sum_rel_abund)), 
                         .before = 0)) %>%
  ungroup() %>%
  select(-sum_rel_abund)

arc.plot = ggplot(arc, aes(x=sample, y=perc_sum_rel_abund*100, fill=phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage relative abundance (reads)",
       x = "Sample",
       fill="Phylum") +
  scale_fill_manual(values = clrs2) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust = 0),
        legend.position = "bottom")



### Figure 3 ###
### alpha diversity ###
sample_by_genus = pivot_wider(reads_by_genus %>% select(sample, genus, rel_abund), names_from = genus, values_from = rel_abund)
sample_by_genus[is.na(sample_by_genus)] = 0

# get shannon diversity
shannon = diversity(select(sample_by_genus, -sample), index = "shannon")
shannon = as.data.frame(shannon)

# get metadata names to look nice
shannon$SampleID = sample_by_genus$sample
shannon = merge(shannon, metadata, on="SampleID")
shannon$Season = str_replace_all(shannon$Season, c("^monsoon$" = "Monsoon\nJul-Aug",
                                                   "^post monsoon$" = "Post-monsoon\nOct-Nov",
                                                   "^winter$" = "Winter\nJan-Feb",
                                                   "^pre monsoon$" = "Pre-monsoon\nApr-May"))
  
shannon$Season = factor(shannon$Season, levels = c("Winter\nJan-Feb","Pre-monsoon\nApr-May","Monsoon\nJul-Aug","Post-monsoon\nOct-Nov"))

farmer_AD = ggplot(shannon, aes(x=Farmer, y=shannon, fill = Farmer)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, fill = "black") +
  scale_fill_manual(values = unique(clrs)) +
  theme_minimal() +
  theme(legend.position="none") +
  labs(y="Shannon Diversity", x = "Farm")

main_crop_AD = ggplot(shannon, aes(x=Main_crop, y=shannon)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, (aes(fill = SampleID))) +
  fillScale +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  labs(x="Crop type")

upazila_AD = ggplot(shannon, aes(x=Upazila, y=shannon, fill=Upazila)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, (aes(fill = SampleID))) +
  fillScale +
  theme_minimal() +
  theme(legend.position="none") +
  labs(y="Shannon Diversity")

season_AD = ggplot(shannon, aes(x=Season, y=shannon)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, (aes(fill = SampleID))) +
  fillScale +
  #geom_text(aes(label = SampleID)) +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

plot_grid( 
  farmer_AD, 
  main_crop_AD, 
  upazila_AD, 
  season_AD,
  labels = c("A","B","C","D"), label_size = 12, 
  align = "h")

### compare other alpha diversity metrics 
shannon = diversity(select(sample_by_genus, -sample), index = "shannon")
simpson = diversity(select(sample_by_genus, -sample), index = "simpson")
inv.simpson = diversity(select(sample_by_genus, -sample), index = "invsimpson")
# chao1 requires whole numbers to multiple by a large number to maintain proportions and round to nearest whole number
chao1 = estimateR(round(select(sample_by_genus, - sample)*1000000000)) 
chao1 = chao1[2,]
evenness = diversity(select(sample_by_genus, -sample)) / 
  log(specnumber(select(sample_by_genus, -sample)))
alpha.div = data.frame(shannon, simpson, inv.simpson, chao1, evenness)
alpha.div$SampleID = sample_by_genus$sample
alpha.div = merge(alpha.div, metadata, on="SampleID")

### run linear (mixed) models to check for significants, Output at Supplimentary table 1 ###
Farmer = summ(lm(shannon ~ Farmer, data=alpha.div), confint = T, digits = 3)
Main_crop = summ(lmer(shannon ~ Main_crop + (1|Farmer), data=alpha.div), confint = T)
Upazila = summ(lmer(shannon ~ Upazila + (1|Farmer), data=alpha.div), confint = T)
Village = summ(lmer(shannon ~ Village + (1|Farmer) + (1|Upazila), data=alpha.div), confint = T)
Season = summ(lmer(shannon ~ Season + (1|Farmer), data=alpha.div), confint = T)




### Figure 4 ###
########## Betadiversity ###########
sample_by_genus = pivot_wider(reads_by_genus %>% 
                                select(sample, genus, rel_abund), 
                              names_from = genus, values_from = rel_abund) %>%
  arrange(sample) %>%
  column_to_rownames(var = "sample")
  
metadata = metadata %>% arrange(SampleID)
sample_by_genus[is.na(sample_by_genus)] = 0

dist_matrix = vegdist(sample_by_genus, method="bray")

mds = metaMDS(dist_matrix, parallel = 8, trymax=100)
sample_by_genus$MDS1 = mds$points[,1]
sample_by_genus$MDS2 = mds$points[,2]

ef = envfit(mds, metadata, permu = 999)
env.factors = rbind(scores(ef, display="vectors"), 
                 scores(ef, display = "factors"))

hull = sample_by_genus %>%
  mutate(Farm = str_sub(rownames(sample_by_genus),2,2)) %>%
  group_by(Farm) %>%
  slice(chull(MDS1,MDS2))

### Get betadiversity by farm ###
permanova = adonis(dist_matrix ~ str_sub(rownames(sample_by_genus),2,2), 
                   method="bray",perm=999)

f.env.var = as.data.frame(env.factors[(str_detect(rownames(env.factors), "Farmer")),])
row.names(f.env.var) = gsub("Farmer","",row.names(f.env.var))
f.env.var$clr = unique(clrs)[1:6]

BD_farm = ggplot(sample_by_genus) +
  geom_point(aes(x = MDS1, y = MDS2,
                 colour = rownames(sample_by_genus)), 
             size=3) +
  colScale +
  geom_polygon(data=hull,alpha = 0.5, aes(x=MDS1,y=MDS2, fill = Farm)) +
  fillScale +
  #coord_fixed() +
  geom_segment(data = f.env.var,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = f.env.var$clr) +
  geom_label(data = f.env.var, 
                   aes(x = NMDS1, y = NMDS2, label = rownames(f.env.var)),
                   size = 4, 
                   fill = f.env.var$clr,
                   colour = c(rep("white", nrow(f.env.var)))) +
  annotate(geom="label", 
           label = paste("PERMANOVA, R2 = ",
                         round(permanova$aov.tab$R2[1], 3),
                         ", p = ", permanova$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 1) + 
  annotate(geom="label", 
           label = paste("envfit R2 = ",round(ef$factors$r["Farmer"], 3),
                         ", p = ", ef$factors$pvals["Farmer"], sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 2) +
  theme_minimal() +
  theme(legend.position="none")

### Get betadiversity by Upazila ###
u.permanova = adonis(dist_matrix ~ metadata$Upazila, method="bray",perm=999)
pairwise.adonis(dist_matrix, metadata$Upazila, p.adjust.m = "BH")

u.env.var = as.data.frame(env.factors[(str_detect(rownames(env.factors), "Upazila")),])
u.env.var$clr = 0
u.env.var["UpazilaJamalpur Sadar",]$clr = clrs$`Jamalpur Sadar`
u.env.var["UpazilaMuktagacha",]$clr = clrs$Muktagacha
u.env.var["UpazilaTarakanda",]$clr = clrs$Tarakanda
row.names(u.env.var) = gsub("Upazila","",row.names(u.env.var))

hull = sample_by_genus %>%
  select(MDS1, MDS2) %>%
  mutate(SampleID = rownames(.)) %>%
  inner_join(metadata, by = "SampleID") %>%
  group_by(Upazila) %>%
  slice(chull(MDS1,MDS2))

BD_upazila = ggplot(sample_by_genus) +
  geom_point(aes(x = MDS1, y = MDS2, 
                 colour = rownames(sample_by_genus)), size=3) +
  colScale +
  geom_polygon(data=hull,alpha = 0.5, aes(x=MDS1,y=MDS2, fill = Upazila)) +
  fillScale +
  #coord_fixed() +
  geom_segment(data = u.env.var,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = u.env.var$clr) +
  geom_label(data = u.env.var, 
                   aes(x = NMDS1, y = NMDS2, label = rownames(u.env.var)),
                   size = 4, 
                   fill = u.env.var$clr,
                   colour = c(rep("white",3))) +
  annotate(geom="label", 
           label = paste("PERMANOVA, R2 = ",
                         round(u.permanova$aov.tab$R2[1], 3),
                         ", p = ", u.permanova$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 1) +
  annotate(geom="label", 
           label = paste("envfit R2 = ",round(ef$factors$r["Upazila"], 3),
                         ", p = ", ef$factors$pvals["Upazila"], sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 2) +
  theme_minimal() +
  theme(legend.position="none")

plot_grid(BD_farm, BD_upazila,
          labels = c("A","B"), 
          label_size = 12,
          ncol = 1)

### run a pairwise PERMANONVA to determine which Upazila is signifcantly different from the other
### from here https://github.com/pmartinezarbizu/pairwiseAdonis
source("C:/Users/agb214/R files/ARG Project/scripts/pairwise.adonis.R")
summary(pairwise.adonis(dist_matrix, metadata$Upazila, p.adjust.m='BH'))

### Figure 5 ###
### KEGG pathway map ###

### removal of human or unknown pathways not in microorganisms
KEGGs = KEGGs[KEGGs$Cat1_Name != "Human Diseases", ]
KEGGs = KEGGs[KEGGs$Cat1_Name != "Not Included in Pathway or Brite", ]
KEGGs = KEGGs[KEGGs$Cat1_Name != "Brite Hierarchies", ]
KEGGs = KEGGs[KEGGs$Cat1_Name != "Organismal Systems", ]
KEGGs = KEGGs[KEGGs$Cat1_Name != "NA", ]

# groupby and summarise by sample and category 
Cat2 <- df %>% 
  dplyr::group_by(SampleID, Cat2_Name) %>%
  summarise(TPM_sum = sum(TPM), .groups = "drop") %>%
  pivot_wider(names_from = Cat2_Name, values_from = TPM_sum) %>%
  drop_na(SampleID) %>%
  column_to_rownames(var = "SampleID") %>%
  replace(is.na(.), 0) %>%
  scale(center = T, scale = T)

Cat2 <- Cat2[, which(colSums(Cat2) != 0)]
Cat2 = Cat2[, order(colSums(-Cat2))]

pheatmap(Cat2, scale = "none", cluster_rows = F)
adonis(Cat2 ~ substr(rownames(Cat2), 1, 2))

### Figure 6 ###
########## ARGs ##########
### setup ### 

clrs = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
         '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', 
         '#800000', '#aaffc3', '#808000', '#000075', '#a9a9a9', '#000000')
names(clrs) = unique(ARGs.expand$Drug.Class)
colScale <- scale_fill_manual(name = "Drug.Class",values = clrs)
######
ARGs_ex.adef = 
  subset(ARGs.expand, Best_Hit_ARO != "adeF")

# plot heatmap of ARGs by Drug.class
metadata_trim = metadata %>%
  select(SampleID, Season, Upazila) %>%
  column_to_rownames(var = "SampleID")

metadata_trim$Season = str_replace_all(metadata_trim$Season, c("^monsoon$" = "Monsoon (Jul-Aug)",
                                                   "^post monsoon$" = "Post-monsoon (Oct-Nov)",
                                                   "^winter$" = "Winter (Jan-Feb)",
                                                   "^pre monsoon$" = "Pre-monsoon (Apr-May)"))


annoCol = list(Upazila = c(`Jamalpur Sadar` = "#e377c2", 
                         Muktagacha = "#17becf",
                         Tarakanda = "#bcbd22"),
               Season = c(`Monsoon (Jul-Aug)` = "#fabed4",
               `Post-monsoon (Oct-Nov)` = "#469990", 
               `Winter (Jan-Feb)` = "#a9a9a9",
               `Pre-monsoon (Apr-May)` = "#000000"))

scaled = ARGs_ex.adef %>% 
  inner_join(metadata, by = "SampleID") %>%
  group_by(Drug.Class, sample, Best_Hit_ARO) %>% 
  summarise(sum = sum(TPM_abundance), .groups = "drop_last") %>%
  summarise(mean = mean(sum), .groups = "drop") %>%
  pivot_wider(names_from = Drug.Class, values_from = mean) %>%
  arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  as.matrix(.) %>%
  replace(is.na(.), 0) %>%
  scale(., center = T, scale = T) %>%
  pheatmap(., scale = "none", 
           annotation_row = metadata_trim,
           annotation_colors = annoCol, 
           cluster_rows = F)

# stats test for ARG enrichement in ponds
### Are specific ponds enriched in ARGs ###
df = ARGs_ex.adef %>% 
  group_by(SampleID, Drug.Class) %>% 
  summarise(sum = sum(TPM_abundance), .groups = "drop") %>% 
  pivot_wider(names_from = Drug.Class, values_from = sum) %>%  
  replace(is.na(.), 0) %>% 
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

with.meta = df %>%
  inner_join(metadata, by = "SampleID") %>%
  arrange(SampleID)

adonis(df ~ with.meta$Farmer, method="bray",perm=999)
adonis(df ~ with.meta$Upazila, method="bray",perm=999)
adonis(df ~ with.meta$Village, method="bray",perm=999)
adonis(df ~ with.meta$Main_crop, method="bray",perm=999)
adonis(df ~ with.meta$Season, method="bray",perm=999)


#### Figure 7 
### Pathogen search ###

pathogen = tibble(read.table("C:/Users/agb214/R files/ARG Project/metadata/ABSA_tabular.csv", sep = ",", header = T))

pathogen_trm = 
  pathogen %>%
  separate(Name, c("genus", "species"), sep = " ", extra = "drop") %>%
  filter(AP == "y")

ARGs_ARO_safe$Best_Hit_ARO = str_replace_all(ARGs_ARO_safe$Best_Hit_ARO, c("Mycobacterium tuberculosis folC with mutation conferring resistance to para-aminosalicylic acid" = "folC", 
                                                                           "Mycobacterium tuberculosis intrinsic murA conferring resistance to fosfomycin" = "murA",
                                                                           "Mycobacterium tuberculosis rpsL mutations conferring resistance to Streptomycin" = "rpsL",
                                                                           "Acinetobacter baumannii AbaQ" = "AbaQ"))


df = filter(ARGs_ARO_safe, grepl(paste(unique(pathogen_trm$genus),"\\b", 
                              sep = "", collapse = "|"), 
                              ARGs_ARO_safe$genus, ignore.case = T)) %>%
  filter(Best_Hit_ARO != "adeF")

ggplot(df, aes(x=sample, y=TPM_abundance, fill=Drug.Class)) +
  geom_col(position = "dodge") +
  facet_wrap(Best_Hit_ARO~species, scales = "free_x",
             labeller = label_wrap_gen(width=25)) +
  colScale +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.4), 
        legend.position = "bottom") +
  guides(fill=guide_legend(title="Conferred resistance")) +
  ylab("Abundance (TPM)") +
  xlab("Sample")

### Supplimentary Figure 3 ###
### ARG abundance by upazila as boxplots ###
ARGs_ex.adef = 
  subset(ARGs.expand, Best_Hit_ARO != "adeF")

ARGs_ex.adef %>% 
  inner_join(metadata, by = "SampleID") %>%
  group_by(Drug.Class, sample, Upazila) %>% 
  summarise(sum = sum(TPM_abundance), .groups = "drop") %>%
  ggplot(., aes(fill = Drug.Class, y = sum, x = Upazila)) +
  #geom_bar(position = "dodge", stat = "identity") +
  geom_boxplot() +
  colScale +
  facet_wrap(~Drug.Class, ncol = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Upazila", y = "Abundance (TPM)")

### Supplimentary Figure 4 ###
### ARG abundance by gene and pond sample ###
ARGs_ex.adef = 
  subset(ARGs.expand, Best_Hit_ARO != "adeF")

ggplot(ARGs_ex.adef, aes(fill=Drug.Class, y=TPM_abundance, x=sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  facet_wrap(~Drug.Class, ncol = 3) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90)) +
  colScale +
  ylab("Abundance (TPM) \nARGs within ponds summarised by antibiotic class") +
  xlab("Sample")

### Supplimentary figure 5 ### 
### plasmid annotation ###

clrs3 = c("#911eb4","#e6194B","#ffe119","#4363d8","#000000","#f58231","#FFFFFF")
names(clrs3) = unique(plas_strict$func)
fillScale <- scale_fill_manual(name = "ID",values = clrs3)


ggplot(plas_strict, aes(xmin=start, xmax=end, y=sample, 
                         fill=ID, label=func, forward = direction)) +
  geom_gene_arrow(arrowhead_height = unit(10, "mm"), 
                  arrowhead_width = unit(5, "mm"),
                  arrow_body_height = unit(10, "mm")) +
  geom_gene_label(align = "left", min.size = 4, height = 10, 
                  padding.x = unit(3, "mm"), colour = "black") +
  facet_wrap(~ sample, scales = "free", ncol = 1) +
  fillScale +
  theme_genes()
# edited in affinity desinger for clarity

### Supplimentary Figure 6 ####
### ARGs by phyla ####
genus_ARGS = 
  ARGs.expand %>% 
  filter(!grepl("Unclassified|unclassified Bacteria phylum", phylum)) %>%
  filter(Best_Hit_ARO != "adeF") %>%
  select(Drug.Class, genus, TPM_abundance) %>%
  group_by(Drug.Class, genus) %>%
  summarise(mean_TPM = mean(TPM_abundance), .groups = "drop_last") %>%
  left_join(., ARGs.expand, by=c("Drug.Class","genus")) %>%
  ungroup %>%
  group_by(phylum, Drug.Class) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

new_clrs = clrs[names(clrs) %in% unique(genus_ARGS$Drug.Class)]
colScale <- scale_fill_manual(name = "Conferred resistance",values = new_clrs)
  
ggplot(genus_ARGS, aes(x=phylum, y=perc*100, fill=Drug.Class)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Percentage of ARGs conferring resistance to different antibiotic drug classes") +
  xlab("Bacterial Phylum") +
  colScale

### Supplimentary Table 2 ###
ARGs_ex.adef %>%
  group_by(SampleID, Drug.Class) %>%
  summarise(Total_TPM_abundance = sum(TPM_abundance), .groups = "drop") %>%
  pivot_wider(names_from = Drug.Class, values_from = Total_TPM_abundance) %>%
  write.csv(., file = "Total ARG transcripts per kilobase million (TPM) abundance summarise by conferred antibiotic resistance class.csv",
            row.names = F)

### Supplimentary Table 2 from 06_recA_abundance.sh line 4

### Supplimentary Table 4 ###
ARGs_ARO_safe %>%
  group_by(Best_Hit_ARO, sample) %>%
  summarise(sum = sum(TPM_abundance), .groups = "drop") %>%
  pivot_wider(names_from = "Best_Hit_ARO", values_from = "sum") %>%
  arrange(sample) %>%
  write_csv(., "ARG_abundance_by_sample.csv")
