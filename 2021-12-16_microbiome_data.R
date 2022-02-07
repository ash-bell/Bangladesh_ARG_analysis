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

set.seed(1234)
setwd("C:/Users/agb214/R files/ARG Project/")

source("scripts/Ponds_by_ARGs_v2.R")
source("scripts/Ponds_by_contigs_v2.R")
source("scripts/Ponds_by_plasmids_v2.R")

### alpha diversity ###
##should this be by count not rel abund???
pond_by_genus = pivot_wider(ctg_rel_abund_genus %>% select(-superkingdom), names_from = genus, values_from = rel_abund)
pond_by_genus = as.data.frame(pond_by_genus)
row.names(pond_by_genus) = pond_by_genus$sample
pond_by_genus$sample = NULL
pond_by_genus[is.na(pond_by_genus)] = 0
#pond_by_genus[pond_by_genus == 0] = NA

alpha_div = tibble(pond = substring(rownames(pond_by_genus),2,2))

shannon = diversity(pond_by_genus, index = "shannon")
shannon = as.data.frame(shannon)
alpha_div$shannon = shannon$shannon
shannon$pond = substring(rownames(shannon),2,2)
shannon.plot = ggplot(shannon, aes(x=pond, y=shannon, fill = pond)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, fill = "black") +
  theme_minimal() +
  theme(legend.position="none")

simpson = diversity(pond_by_genus, index = "simpson")
simpson = as.data.frame(simpson)
alpha_div$simpson = simpson$simpson
simpson$pond = substring(rownames(simpson),2,2)
simpson.plot = ggplot(simpson, aes(x=pond, y=simpson, fill = pond)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, fill = "black") +
  theme_minimal() +
  theme(legend.position="none")

inv.simpson = diversity(pond_by_genus, index = "invsimpson")
inv.simpson = as.data.frame(inv.simpson)
alpha_div$inv.simpson = inv.simpson$inv.simpson
inv.simpson$pond = substring(rownames(inv.simpson),2,2)
inv.simpson.plot = ggplot(inv.simpson, aes(x=pond, y=inv.simpson, fill = pond)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, fill = "black") +
  theme_minimal() +
  theme(legend.position="none")

pielous.evenness = diversity(pond_by_genus, index = "shannon") / log(specnumber(pond_by_genus))
pielous.evenness = as.data.frame(pielous.evenness)
alpha_div$pielous.evenness = pielous.evenness$pielous.evenness
pielous.evenness$pond = substring(rownames(pielous.evenness),2,2)
pielous.evenness.plot = ggplot(pielous.evenness, aes(x=pond, y=pielous.evenness, fill = pond)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, fill = "black") +
  theme_minimal() +
  theme(legend.position="none")

plot_grid(shannon.plot, simpson.plot, inv.simpson.plot, pielous.evenness.plot,
          labels = c("A","B","C","D"), label_size = 12)

# Find if any pond alpha diversity is significantly different
# Anova (overall) and Tukey's Honest Significant Difference (HSD) (individual)
options(scipen=999)



alpha_div %>%
  mutate(pond = as_factor(pond)) %>%
  select(pond, shannon) %>%
  leveneTest(shannon ~ pond, data = .,)

alpha_div %>% 
  mutate(pond = as_factor(pond)) %>%
  pivot_longer(!pond, names_to = "div_metric", values_to = "alpha_div") %>%
  group_by(div_metric) %>%
#  group_map(~ leveneTest(.x$alpha_div ~ pond, data = .)) %>%
  group_map(~ shapiro.test(x=residuals(object=aov(.x$alpha_div ~ pond, data= .))))

welch.res.shannon = oneway.test(shannon ~ pond, data = alpha_div)
aov.res.simpson = aov(simpson ~ pond, data = alpha_div)
aov.res.inv.simpson = aov(inv.simpson ~ pond, data = alpha_div)
kruskal.res.pielous = kruskal.test(pielous.evenness ~ pond, data = alpha_div)


summary(aov.res.simpson)
summary(aov.res.inv.simpson)

GH.shannon = games_howell_test(shannon ~ pond, data = alpha_div)
Tukey.shannon = TukeyHSD(welch.res.shannon)
Tukey.inv.simpson = TukeyHSD(aov.res.inv.simpson)
WC.pielous = pairwise.wilcox.test(alpha_div$pielous.evenness, alpha_div$pond, p.adjust.method = "BH")

any(GH.shannon$p.adj < 0.05)
any(Tukey.simpson$pond[,4] < 0.05)
any(Tukey.inv.simpson$pond[,4] < 0.05)
any(WC.pielous$p.value < 0.05)

### rarefraction ###
pond_by_genus = pivot_wider(genus_count %>% select(-superkingdom), names_from = genus, values_from = count)
pond_by_genus = as.data.frame(pond_by_genus)
row.names(pond_by_genus) = pond_by_genus$sample
pond_by_genus$sample = NULL
pond_by_genus[is.na(pond_by_genus)] = 0
#pond_by_genus[pond_by_genus == 0] = NA



#Srar <- rarefy(round(pond_by_genus*1000000), min(rowSums(round(pond_by_genus*1000000))))
#rarecurve(round(pond_by_genus[, 1:1000]*10), sample=min(rowSums(round(pond_by_genus[, 1:1000]*10))))

########## abundance barplot #####

class_abundance = contigs %>% 
  filter(!genus == "unclassified") %>%
  summarise(sum = sum(TPM_abundance))

unclass_abundance = contigs %>% 
  filter(genus == "unclassified") %>%
  summarise(sum = sum(TPM_abundance))

total_abund = contigs %>% 
  summarise(sum = sum(TPM_abundance))

ctg.rel.abund = contigs %>% 
  filter(!genus == "unclassified") %>%
  group_by(sample, genus) %>%
  summarise(mean_TPM = mean(TPM_abundance), .groups = "drop_last") %>%
  left_join(., contigs, by = c("sample","genus")) %>%
  distinct(sample,genus,mean_TPM, .keep_all = TRUE) %>%
  bind_rows(., filter(.data = contigs, genus == "unclassified"))

ctg.rel.abund[ctg.rel.abund$genus == "unclassified", ]$mean_TPM = 
  ctg.rel.abund[ctg.rel.abund$genus == "unclassified", ]$TPM_abundance
  
ctg.rel.abund = ctg.rel.abund %>% group_by(sample) %>% 
  mutate(rel_abund = mean_TPM/sum(mean_TPM)) %>%
  select(-TPM_abundance, -mean_TPM) %>%
  ungroup()

bac = ctg.rel.abund %>% 
  filter(superkingdom == "Bacteria") %>%
  select(sample,phylum,rel_abund) %>%
  group_by(sample, phylum) %>%
  summarise(sum_rel_abund = sum(rel_abund), .groups = "drop_last") %>%
  mutate(perc_sum_rel_abund = sum_rel_abund/sum(sum_rel_abund)) %>%
  slice_max(sum_rel_abund, n=10) %>%
  group_modify(~ add_row(.x, phylum = "Other", 
                         perc_sum_rel_abund = (1-sum(.$perc_sum_rel_abund)), 
                         .before = 0)) %>%
  ungroup() %>%
  select(-sum_rel_abund)

bac.plot = ggplot(bac, aes(x=sample, y=perc_sum_rel_abund, fill=phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage relative abundance (TPM)",
       x = "Sample",
       fill="Phylum")

vir = ctg.rel.abund %>% 
  filter(superkingdom == "Viruses") %>%
  select(sample,phylum,rel_abund) %>%
  group_by(sample, phylum) %>%
  summarise(sum_rel_abund = sum(rel_abund), .groups = "drop_last") %>%
  mutate(perc_sum_rel_abund = sum_rel_abund/sum(sum_rel_abund)) %>%
  slice_max(sum_rel_abund, n=10) %>%
  group_modify(~ add_row(.x, phylum = "Other", 
                         perc_sum_rel_abund = (1-sum(.$perc_sum_rel_abund)), 
                         .before = 0)) %>%
  ungroup() %>%
  select(-sum_rel_abund)

vir.plot = ggplot(vir, aes(x=sample, y=perc_sum_rel_abund, fill=phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage relative abundance (TPM)",
       x = "Sample",
       fill="Phylum")

arc = ctg.rel.abund %>% 
  filter(superkingdom == "Archaea") %>%
  select(sample,phylum,rel_abund) %>%
  group_by(sample, phylum) %>%
  summarise(sum_rel_abund = sum(rel_abund), .groups = "drop_last") %>%
  mutate(perc_sum_rel_abund = sum_rel_abund/sum(sum_rel_abund)) %>%
  slice_max(sum_rel_abund, n=10) %>%
  group_modify(~ add_row(.x, phylum = "Other", 
                         perc_sum_rel_abund = (1-sum(.$perc_sum_rel_abund)), 
                         .before = 0)) %>%
  ungroup() %>%
  select(-sum_rel_abund)

arc.plot = ggplot(arc, aes(x=sample, y=perc_sum_rel_abund, fill=phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage relative abundance (TPM)",
       x = "Sample",
       fill="Phylum")

plot_grid(bac.plot, arc.plot, vir.plot,
          labels = c("A","B","C"), label_size = 12,
          ncol = 1)


df %>% group_by(Drug.Class, sample) %>%
  summarise(TPM_abundance = sum(TPM_abundance), .groups = "drop_last") %>%
  inner_join(., metadata, by=c("sample" = "SampleID")) %>%
  select(Drug.Class, TPM_abundance, Main_crop) %>%
  #group_map(~ leveneTest(TPM_abundance ~ Site, data = .))
  #group_map(~ shapiro.test(x=residuals(object=aov(TPM_abundance ~ Drug.Class, data= .))))
  #group_map(~ kruskal.poss(TPM_abundance ~ Site, data = .))
  #group_map(~ kruskal.poss(TPM_abundance ~ Main_crop, data = .))
  group_map(~ wilcox.poss(.$TPM_abundance, .$Main_crop, p.adjust.method = "BH"))

kruskal.poss = possibly(.f = kruskal.test, otherwise = "Error")
wilcox.poss = possibly(.f = pairwise.wilcox.test, otherwise = "Error")

vir %>% mutate(farm = substr(sample,2,2)) %>%
  mutate(district = ifelse(farm == c(1,2,3), "A", 
                       ifelse(farm == c(4,8), "B", "C"))) %>%
  select(-sample, farm) %>%
  group_by(phylum) %>% 
  #group_map(~ kruskal.poss(perc_sum_rel_abund ~ district, data = .))
  group_map(~ wilcox.poss(.$perc_sum_rel_abund, .$district, p.adjust.method = "BH"))

  summarise(abund = sum(perc_sum_rel_abund), .groups = "drop") %>% mutate(abund = abund/sum(abund)*100)

########## Betadiversity ###########
ctg.rel.abund = contigs %>% 
  filter(!genus == "unclassified") %>%
  group_by(sample, genus) %>%
  summarise(mean_TPM = mean(TPM_abundance), .groups = "drop_last") %>%
  left_join(., contigs, by = c("sample","genus")) %>%
  distinct(sample,genus,mean_TPM, .keep_all = TRUE) %>%
  bind_rows(., filter(.data = contigs, genus == "unclassified"))

ctg.rel.abund[ctg.rel.abund$genus == "unclassified", ]$mean_TPM = 
  ctg.rel.abund[ctg.rel.abund$genus == "unclassified", ]$TPM_abundance

ctg.rel.abund = ctg.rel.abund %>% group_by(sample) %>% 
  mutate(rel_abund = mean_TPM/sum(mean_TPM)) %>%
  select(-TPM_abundance, -mean_TPM) %>%
  ungroup()

df = ctg.rel.abund %>% 
  select(-Contig,-superkingdom,-phylum,-class,-order,-family,-species,-pond) %>%
  filter(!genus == "unclassified") %>%
  pivot_wider(names_from = genus, values_from = rel_abund) %>%
  replace(is.na(.), 0) %>%
  remove_rownames(.) %>% 
  column_to_rownames(., var = "sample")

dist_matrix = vegdist(df, method="bray")

mds = metaMDS(dist_matrix, parallel = 4, k=2, trymax=100)
df$MDS1 = mds$points[,1]
df$MDS2 = mds$points[,2]

metadata = read.csv("../metadata/BD_fish_ponds_pooled_samples_for_metagenome_seq.csv")

ef = envfit(mds, metadata, permu = 999)
env.factors = rbind(scores(ef, display="vectors"), 
                 scores(ef, display = "factors"))

hull = df %>%
  mutate(Site = str_sub(rownames(df),2,2)) %>%
  group_by(Site) %>%
  slice(chull(MDS1,MDS2))

### by pond ###
clrs = c(rep('#e6194B', 4),
         rep('#3cb44b', 4),
         rep('#ffe119', 4),
         rep('#4363d8', 4),
         rep('#f58231', 4),
         rep('#911eb4', 4))
names(clrs) = str_sub(rownames(df),2,2)
colScale <- scale_colour_manual(name = "Site",values = clrs)
fillScale <- scale_fill_manual(name = "Site",values = clrs)

permanova = adonis(dist_matrix ~ str_sub(rownames(df),2,2), method="bray",perm=999)
permutation_test = permutest(betadisper(dist_matrix , str_sub(rownames(df),2,2)), pairwise = T)

env.var = as.data.frame(env.factors[(str_detect(rownames(env.factors), "Farmer")),])
env.var$clr = unique(clrs)

ggplot(df) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = str_sub(rownames(df),2,2)), size=3) +
  colScale +
  geom_polygon(data=hull,alpha = 0.5, aes(x=MDS1,y=MDS2, fill = Site)) +
  fillScale +
  coord_fixed() +
  geom_segment(data = env.var,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = env.var$clr) +
  geom_label_repel(data = env.var, 
                   aes(x = NMDS1, y = NMDS2, label = rownames(env.var)),
                   size = 3, 
                   fill = env.var$clr,
                   colour = c(rep("white", nrow(env.var)))) +
  annotate(geom="label", 
           label = paste("PERMANOVA, R^2 = ",
                         round(permanova$aov.tab$R2[1], 3),
                         ", p = ", permanova$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 1) +
  annotate(geom="label", 
           label = paste("PERMDISP2, p = ",
                         permutation_test$tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 2)

### by Upazila ###
u.permanova = adonis(dist_matrix ~ metadata$Upazila, method="bray",perm=999)
u.permutation_test = permutest(betadisper(dist_matrix , metadata$Upazila), pairwise = T)

env.var = as.data.frame(env.factors[(str_detect(rownames(env.factors), "Upazila")),])
env.var$clr = c("red","blue","darkgreen")

clrs = c(rep('darkblue', 4),
         rep('deepskyblue', 4),
         rep('aquamarine', 4),
         rep('red', 4),
         rep('limegreen', 4),
         rep('darkred', 4))
names(clrs) = str_sub(rownames(df),2,2)
colScale <- scale_colour_manual(name = "Site",values = clrs)
fillScale <- scale_fill_manual(name = "Site",values = clrs)

ggplot(df) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = str_sub(rownames(df),2,2)), size=3) +
  colScale +
  geom_polygon(data=hull,alpha = 0.5, aes(x=MDS1,y=MDS2, fill = Site)) +
  fillScale +
  coord_fixed() +
  geom_segment(data = env.var,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = env.var$clr) +
  geom_label_repel(data = env.var, 
                   aes(x = NMDS1, y = NMDS2, label = rownames(env.var)),
                   size = 3, 
                   fill = env.var$clr,
                   colour = c(rep("white",3))) +
  annotate(geom="label", 
           label = paste("PERMANOVA, R^2 = ",
                         round(u.permanova$aov.tab$R2[1], 3),
                         ", p = ", u.permanova$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 1) +
  annotate(geom="label", 
           label = paste("PERMDISP2, p = ",
                         u.permutation_test$tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 2) +
  annotate(geom="label", 
           label = paste("Goodness of fit R^2 = 0.2313, p = 0.039"), 
           x=Inf, y=Inf, hjust = 1, vjust = 3)


########## ARGs ##########
### setup ### 

clrs = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
         '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', 
         '#800000', '#aaffc3', '#808000', '#000075', '#a9a9a9', '#000000')
names(clrs) = unique(ARGs.expand$Drug.Class)
colScale <- scale_fill_manual(name = "Drug.Class",values = clrs)

########## Plasmid identity, abundance, location, resistance ##########
ggplot(plas_strict, aes(x=sample, y=TPM_abundance, fill=Drug.Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~species, scales = "free_x", 
             labeller = label_wrap_gen(width=30)) +
  ylab("Abundance (TPM)") +
  colScale +
  guides(fill=guide_legend(title="Conferred resistance"))

################ Bacterial #################
#### without adeF gene ###
ARGs_ex.adef = 
  subset(ARGs.expand, Best_Hit_ARO != "adeF")

ggplot(ARGs_ex.adef, aes(fill=Drug.Class, y=TPM_abundance, x=sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  facet_wrap(~Drug.Class) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=90)) +
  colScale +
  ylab("Abundance (TPM) \nARGs within ponds summarised by antibiotic class") +
  xlab("Ponds")

ARGs_ARO_safe$Best_Hit_ARO = str_replace_all(ARGs_ARO_safe$Best_Hit_ARO, c("Mycobacterium tuberculosis folC with mutation conferring resistance to para-aminosalicylic acid" = "folC", 
                                                                       "Mycobacterium tuberculosis intrinsic murA conferring resistance to fosfomycin" = "murA",
                                                                       "Mycobacterium tuberculosis rpsL mutations conferring resistance to Streptomycin" = "rpsL",
                                                                       "Acinetobacter baumannii AbaQ" = "AbaQ"))

ARGs_ARO_safe = subset(ARGs_ARO_safe, Best_Hit_ARO != "adeF")

print(ARGs_ARO_safe %>% 
  group_by(Best_Hit_ARO, sample) %>% 
  summarise(n = n()) %>% 
  summarise(n = n()), n=100)

ggplot(ARGs_ARO_safe, aes(fill=Drug.Class, y=TPM_abundance, x=sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  facet_wrap(~Best_Hit_ARO) +
  theme(axis.text.x=element_text(angle=90)) +
  colScale +
  ylab("Abundance (TPM) \nARGs within ponds summarised by ARG") +
  xlab("Ponds")


### Are specific ponds enriched in ARGs ###
df = ARGs.expand %>% filter(Best_Hit_ARO != "adeF") %>%
  select(Drug.Class, TPM_abundance, sample) %>%
  mutate(Site = str_sub(sample,2,2)) %>%
  mutate(Upazila = ifelse(Site == 1 | Site == 2 | Site == 3, "Muktagacha (F1,2,3)",
                          ifelse(Site == 4 | Site == 8, "Jamalpur Sadar (F4,8)",
                                 ifelse(Site == 6, "Tarakanda (F6)", NA))))

ggplot(df, aes(x=Drug.Class, y=TPM_abundance, fill = Drug.Class)) +
  geom_boxplot() +
  colScale +
  facet_wrap(~Site) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

metadata = read.csv("../metadata/BD_fish_ponds_pooled_samples_for_metagenome_seq.csv")

kruskal.poss = possibly(.f = kruskal.test, otherwise = "Error")
wilcox.poss = possibly(.f = pairwise.wilcox.test, otherwise = "Error")

df %>% group_by(Drug.Class, sample) %>%
  summarise(TPM_abundance = sum(TPM_abundance), .groups = "drop_last") %>%
  inner_join(., metadata, by=c("sample" = "SampleID")) %>%
  select(Drug.Class, TPM_abundance, Main_crop) %>%
  #group_map(~ leveneTest(TPM_abundance ~ Site, data = .))
  #group_map(~ shapiro.test(x=residuals(object=aov(TPM_abundance ~ Drug.Class, data= .))))
  #group_map(~ kruskal.poss(TPM_abundance ~ Site, data = .))
  #group_map(~ kruskal.poss(TPM_abundance ~ Main_crop, data = .))
  group_map(~ wilcox.poss(.$TPM_abundance, .$Main_crop, p.adjust.method = "BH"))

tmp = df %>% group_by(Drug.Class, sample) %>%
  summarise(TPM_abundance = sum(TPM_abundance), .groups = "drop_last") %>%
  mutate(Site = str_sub(sample,2,2))

welch.res.shannon = oneway.test(shannon ~ pond, data = alpha_div)
aov.res.simpson = aov(simpson ~ pond, data = alpha_div)
aov.res.inv.simpson = aov(inv.simpson ~ pond, data = alpha_div)
kruskal.res.pielous = 


summary(aov.res.simpson)
summary(aov.res.inv.simpson)

GH.shannon = games_howell_test(shannon ~ pond, data = alpha_div)
Tukey.shannon = TukeyHSD(welch.res.shannon)
Tukey.inv.simpson = TukeyHSD(aov.res.inv.simpson)
WC.pielous = pairwise.wilcox.test(alpha_div$pielous.evenness, alpha_div$pond, p.adjust.method = "BH")

any(GH.shannon$p.adj < 0.05)
any(Tukey.simpson$pond[,4] < 0.05)
any(Tukey.inv.simpson$pond[,4] < 0.05)
any(WC.pielous$p.value < 0.05)

### rarefraction ###
pond_by_genus = pivot_wider(genus_count %>% select(-superkingdom), names_from = genus, values_from = count)
pond_by_genus = as.data.frame(pond_by_genus)
row.names(pond_by_genus) = pond_by_genus$sample
pond_by_genus$sample = NULL
pond_by_genus[is.na(pond_by_genus)] = 0
#pond_by_genus[pond_by_genus == 0] = NA








drugclass = ARGs.expand %>%
  select(Drug.Class, sample, TPM_abundance) %>%
  pivot_wider(names_from = Drug.Class, values_from = TPM_abundance, values_fn = sum) %>%
  replace(is.na(.), 0) %>%
  remove_rownames(.) %>% 
  column_to_rownames(., var = "sample")

DC.dis.matrix = vegdist(drugclass, method="bray")

mds = metaMDS(DC.dis.matrix, parallel = 4, k=2, trymax=100)

ef = envfit(mds, drugclass, permu = 999)

drugclass$MDS1 = mds$points[,1]
drugclass$MDS2 = mds$points[,2]

hull = drugclass %>%
  mutate(Site = str_sub(rownames(drugclass),2,2)) %>%
  group_by(Site) %>%
  slice(chull(MDS1,MDS2))

clrs = c(rep('#e6194B', 4),
         rep('#3cb44b', 4),
         rep('#ffe119', 4),
         rep('#4363d8', 4),
         rep('#f58231', 4),
         rep('#911eb4', 4))
names(clrs) = str_sub(rownames(drugclass),2,2)
colScale <- scale_colour_manual(name = "Site",values = clrs)
fillScale <- scale_fill_manual(name = "Site",values = clrs)

env.factors = rbind(scores(ef, display="vectors"), 
                    scores(ef, display = "factors"))

env.var = as.data.frame(env.factors)
env.var$clr = unique(clrs)

ggplot(drugclass) +
  geom_point(mapping=aes(x=MDS1,y=MDS2,colour=str_sub(rownames(drugclass),2,2)),size=3) +
  colScale +
  geom_polygon(data=hull,alpha = 0.5, aes(x=MDS1,y=MDS2, fill=Site)) +
  fillScale +
  coord_fixed() +
  geom_segment(data = env.var,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = env.var$clr) +
  geom_label_repel(data = env.var, 
                   aes(x = NMDS1, y = NMDS2, label = rownames(env.var)),
                   size = 3, 
                   fill = env.var$clr,
                   colour = c(rep("white", nrow(env.var)))) 

+
  annotate(geom="label", 
           label = paste("PERMANOVA, R^2 = ",
                         round(permanova$aov.tab$R2[1], 3),
                         ", p = ", permanova$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 1) +
  annotate(geom="label", 
           label = paste("PERMDISP2, p = ",
                         permutation_test$tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 2)




  




tmp = tmp[!tmp$Drug.Class == "para-aminosalicylic acid", ]
tmp = tmp[!tmp$Drug.Class == "cephamycin", ]

df_aov = tmp %>%
  group_by(Drug.Class) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% 
           map(.x = ., .f = ~ TukeyHSD(aov(TPM_abundance ~ Site, data = .))))

for (i in 1:(length(df_aov$aov_results))) {
  print(df_aov$Drug.Class[i])
  print(df_aov$aov_results[[i]]$Site[,4][df_aov$aov_results[[i]]$Site[,4] < 0.05])
  }

df_aov = tmp %>%
  group_by(Drug.Class) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% 
           map(.x = ., .f = ~ TukeyHSD(aov(TPM_abundance ~ Upazila, data = .))))

for (i in 1:(length(df_aov$aov_results))) {
  print(df_aov$Drug.Class[i])
  print(df_aov$aov_results[[i]]$Upazila[,4][df_aov$aov_results[[i]]$Upazila[,4] < 0.05])
}


##### which taxa responsible ####
genus_ARGS = 
  ARGs.expand %>% 
  select(Drug.Class, genus, TPM_abundance) %>%
  group_by(Drug.Class, genus) %>%
  summarise(mean_TPM = mean(TPM_abundance), .groups = "drop_last") %>%
  left_join(., ARGs.expand, by=c("Drug.Class","genus")) %>%
  ungroup
  
ggplot(genus_ARGS, aes(x=phylum, y=mean_TPM, fill=Drug.Class)) +
  geom_bar(position = "stack", stat = "identity") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
  #facet_grid(~phylum, space = "free_x", scales = "free_x") +
  coord_flip() +
  colScale


sample = slice_sample(ARGs.expand, n=10)
tree = read.tree(text = "((1760988:1,133552:1,2026199:1)1224:1,((2026780:1,2026779:1)203682:1,1032527:1)1783257:1,1504320:1);")
plotTree(tree)

df = ARGs.expand %>% filter(Best_Hit_ARO != "adeF") %>%
  #filter(!superkingdom == "Unclassified") %>%
  select(genus, TPM_abundance, sample) %>%
  pivot_wider(names_from = sample, values_from = TPM_abundance, values_fn = mean) %>%
  replace(is.na(.), 0) %>%
  remove_rownames(.) %>% 
  column_to_rownames(., var = "genus")

dist_matrix = vegdist(df, method="bray")

mds = metaMDS(dist_matrix, parallel = 4, k=3, trymax=100)

metadata = ARGs.expand %>%
  select(Drug.Class, TPM_abundance, genus) %>%
  pivot_wider(names_from = Drug.Class, values_from = TPM_abundance, values_fn = sum) %>%
  replace(is.na(.), 0) %>%
  remove_rownames(.) %>% 
  column_to_rownames(., var = "genus")



##############
### Pathogen search ###

pathogen = tibble(read.table("../metadata/ABSA_tabular.csv", sep = ",", header = T))

pathogen_trm = 
  pathogen %>%
  separate(Name, c("genus", "species"), sep = " ", extra = "drop") %>%
  filter(AP == "y")

df = filter(ARGs_ARO_safe, grepl(paste(unique(pathogen_trm$genus),"\\b", 
                              sep = "", collapse = "|"), 
                              ARGs_ARO_safe$genus, ignore.case = T))

ggplot(df, aes(x=sample, y=TPM_abundance, fill=Drug.Class)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(Best_Hit_ARO~species, scales = "free_x") +
  colScale +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.4))
  #scale_x_discrete(guide = guide_axis(n.dodge=2))




dist_matrix_genus = vegdist(metadata, method="bray")

mds = metaMDS(dist_matrix_genus, parallel = 4, k=3, trymax=100)

ef = envfit(mds, metadata, permu = 999)
ef_genus = envfit(mds, df, permu = 999)

adonis(dist_matrix_genus ~ df$`unclassified Actinomycetia genus`, method="bray",perm=999)
adonis(dist_matrix_genus ~ df$Mycobacterium, method="bray",perm=999)
adonis(dist_matrix_genus ~ df$`unclassified Ilumatobacteraceae genus`, method="bray",perm=999)
adonis(dist_matrix_genus ~ df$`unclassified Mycobacteriaceae genus`, method="bray",perm=999)

env.factors = as.data.frame(scores(ef, display = "vector"))


df$MDS1 = mds$points[,1]
df$MDS2 = mds$points[,2]

clrs = c(rep('#e6194B', 4),
         rep('#3cb44b', 4),
         rep('#ffe119', 4),
         rep('#4363d8', 4),
         rep('#f58231', 4),
         rep('#911eb4', 4))
names(clrs) = str_sub(rownames(df),2,2)
colScale <- scale_colour_manual(name = "Site",values = clrs)

# by Drug.classes
sulfonamide = adonis(dist_matrix ~ metadata$sulfonamide, method="bray",perm=999)
aminoglycoside = adonis(dist_matrix ~ metadata$aminoglycoside, method="bray",perm=999)
acridine = adonis(dist_matrix ~ metadata$`acridine dye`, method="bray",perm=999)

#permutation_test = permutest(betadisper(dist_matrix , metadata$sulfonamide), pairwise = F)

ggplot(df) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = str_sub(rownames(df),2,2)), size=3) +
  colScale +
  coord_fixed() +
  geom_segment(data = env.factors,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), 
               ) +
  geom_text_repel(data = env.factors, 
                  aes(x = NMDS1, y = NMDS2, label = rownames(env.factors)),
                  size = 3, 
                  ) +
  annotate(geom="label", 
           label = paste("Sulfonamide PERMANOVA, R^2 = ",
                         round(sulfonamide$aov.tab$R2[1], 3),
                         ", p = ", sulfonamide$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 1) +
  annotate(geom="label", 
           label = paste("Aminoglycoside PERMANOVA, R^2 = ",
                         round(aminoglycoside$aov.tab$R2[1], 3),
                         ", p = ", aminoglycoside$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 2) +
  annotate(geom="label", 
           label = paste("Acridine dye PERMANOVA, R^2 = ",
                         round(acridine$aov.tab$R2[1], 3),
                         ", p = ", acridine$aov.tab$`Pr(>F)`[1], 
                         sep = ""), 
           x=Inf, y=Inf, hjust = 1, vjust = 3)

