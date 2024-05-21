library(ggplot2)
library(tidyverse)
library(scales)
library(cowplot)
library(patchwork)

meta <- read.table("data/multiome_geo/metadata_csv/h5ad_file_metadata.csv", sep = ",", header = TRUE)

ncells_donor1 <- meta %>% 
  filter((Site == "site1" & DonorNumber == "donor1")) %>% 
  nrow(.)


meta %>% 
  filter(
    (DonorNumber == "donor1" & Site == "site1") | 
      (DonorNumber == "donor2" & Site == "site1") | 
      (DonorNumber == "donor8" & Site == "site4") | 
      (DonorNumber == "donor5" & Site == "site2")
  ) %>% 
  group_by(cell_type) %>% 
  nrow(.)
  summarize(count = n()) %>% 
  View(.)

ncells_donor1_donor2 <- meta %>% 
  filter(
    (DonorNumber == "donor1" & Site == "site1") | 
      (DonorNumber == "donor2" & Site == "site1") | 
      (DonorNumber == "donor8" & Site == "site4") | 
      (DonorNumber == "donor5" & Site == "site2")
  ) %>% 
  group_by(cell_type) %>% 
  summarize(count = n())

p1 <- meta %>% 
  group_by(cell_type) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = reorder(cell_type, -count), y = log10(count), fill=cell_type)) + geom_bar(stat="identity") + scale_y_continuous(labels=scales::label_math()) +
  geom_hline(yintercept = log10(50), color = "blue")  +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank()) +
  labs(y = "Cell barcodes", title = "Cell types, multiome bone marrow, GEO Series GSE194122", subtitle = paste0("All donors, all sites (", format(nrow(meta), big.mark = ","), " cells)"))
  

p2 <- meta %>% 
  # filter((Site == "site1" & DonorNumber == "donor1") | (Site == "site1" & DonorNumber == "donor2") | (Site == "site4" & DonorNumber == "donor8")) %>% 
  filter((Site == "site1" & DonorNumber == "donor1")) %>% 
  group_by(cell_type) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = reorder(cell_type, -count), y = log10(count), fill=cell_type)) + geom_bar(stat="identity") + scale_y_continuous(labels=scales::label_math()) +
  geom_hline(yintercept = log10(50), color = "blue")  +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank()) +
  labs(y = "Cell barcodes", title = "Cell types, multiome bone marrow, GEO Series GSE194122", subtitle = paste0("Donor 1, site 1 (", format(ncells_donor1, big.mark = ","), " cells)"))
  
p2 / p1

p3 <- meta.subset %>% 
  group_by(cell_type) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = reorder(cell_type, -count), y = log10(count), fill=cell_type)) + geom_bar(stat="identity") + scale_y_continuous(labels=scales::label_math()) +
  geom_hline(yintercept = log10(50), color = "blue")  +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank()) +
  labs(y = "Cell barcodes", title = "Cell types, multiome bone marrow, GEO Series GSE194122", subtitle = paste0("~2,000 cells per cell type (", format(nrow(meta.subset), big.mark = ","), " cells)"))

p1 / p2 / p3
