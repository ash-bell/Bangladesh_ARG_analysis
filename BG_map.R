############ World Map #####################
library("tidyverse")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
library("sp")
library(maps)
library(ggrepel)
library(cowplot)

metadata = read.csv("~/../R files/ARG Project/metadata/BD_fish_ponds_pooled_samples_for_metagenome_seq.csv")

world <- ne_countries(scale = "medium", returnclass = "sf")

metadata = 
  metadata %>%
  mutate(Upazila = gsub("Jamalpur_Sadar", "Jamalpur Sadar", metadata$Upazila)) %>%
  mutate(latitude = ifelse(Village == "Kandulia", 24.893648644007897,
                           ifelse(Village == "Malotipur", 24.780456,
                                  ifelse(Village == "Bamunpara", 24.899441,
                                         ifelse(Village == "Gopalpur", 24.561602,
                                                ifelse(Village == "Fulbaria", 24.631393, NA)))))) %>%
  mutate(longitude = ifelse(Village == "Kandulia", 90.76682472089583,
                          ifelse(Village == "Malotipur", 90.284227,
                                 ifelse(Village == "Bamunpara", 90.022605,
                                        ifelse(Village == "Gopalpur", 89.927654,
                                               ifelse(Village == "Fulbaria", 90.269873, NA)))))) %>%
  mutate(Village = ifelse(Village == "Kandulia", "Kandulia\n(F1, F2)",
                           ifelse(Village == "Malotipur", "Malatipur\n(F3)",
                                  ifelse(Village == "Bamunpara", "Bamunpara\n(F4)",
                                         ifelse(Village == "Gopalpur", "Gopalpur\n(F6)",
                                                ifelse(Village == "Fulbaria", "Fulbaria\n(F8)", NA))))))

village = distinct(metadata, Village, .keep_all = T)

metadata = st_as_sf(metadata, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")



#24.780456, 90.284227 Malotipur Malatipur
#24.899441, 90.022605 Bamunpara 
#24.561602, 89.927654 Gopalpur
#24.631393, 90.269873 Fulbaria

#24.764208, 90.256995 Muktagacha 
#24.872342, 90.418575 Tarakanda (Phulpur)
#24.844071, 90.002980 Jamalpur Sadar




Upazila = read_sf("~/../Downloads/BGD_adm/BGD_adm2.shp")

Upazila2 = read_sf("~/../Downloads/administrative_boundary_of_bangladesh_upazila_level/administrative_boundary_of_bangladesh_upazila_level.shp")



ggplot(data=Upazila2) +
  geom_sf(aes(fill = ifelse(NAME_4=="Jamalpur S.", 'Jamalpur S', "green")))


p1 = 
ggplot(data = world) +
  geom_sf() +
  geom_sf(data = Upazila2, aes(fill = ifelse(NAME_4=="Muktagachha", 'Muktagacha',
                                             ifelse(NAME_4=="Jamalpur S.", 'Jamalpur Sadar',
                                                    ifelse(NAME_4=="Phulpur", 'Tarakanda',
                                                           ifelse(NAME_3=="Nasirabad", 'Mymensingh', NA)))))) +
  scale_fill_manual(values=c('#DD4B3E', '#1EA362', '#F6CF65', '#4A89F3'), name="Regions") +
  coord_sf(xlim = c(88, 93), ylim = c(21, 27)) +
  theme(legend.position = "none")




p2 = 
ggplot(data = world) +
  geom_sf() +
  geom_sf(data = Upazila2, aes(fill = ifelse(NAME_4=="Muktagachha", 'Muktagacha',
                                             ifelse(NAME_4=="Jamalpur S.", 'Jamalpur Sadar',
                                                    ifelse(NAME_4=="Phulpur", 'Tarakanda',
                                                           ifelse(NAME_3=="Nasirabad", 'Mymensingh', NA)))))) +
  geom_sf(data = metadata) +
  geom_label_repel(data = village, aes(x=longitude, y=latitude, label=Village), force = 100) +
  scale_fill_manual(values=c('#DD4B3E', '#1EA362', '#F6CF65', '#4A89F3'), name="Regions") +
  coord_sf(xlim = c(89.75, 91.25), ylim = c(24.25,25.25)) +
  theme(legend.position = "bottom")


plot = plot_grid(p1, p2, labels = c("A","B"),
          rel_widths = c(1,3),
          rel_heights = c(3,1))
