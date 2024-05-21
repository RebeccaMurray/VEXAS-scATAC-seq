
# Load the object
obj.multi <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/multiome_geo/multi_ref.rds")

# Select for some donor barcodes
obj.multi$BC <- rownames(obj.multi@meta.data)
obj.multi$to_keep <- grepl("*-(s1d1|s1d2|s4d8|s2d5)$", obj.multi$BC)
table(obj.multi$to_keep)

# Subset
obj.multi.subset <- subset(obj.multi, subset = to_keep)


# Save
saveRDS(obj.multi.subset, "~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/multiome_geo/multi_ref_s1d1_s1d2_s4d8_s2d5.rds")


# V2

# Load the object
obj.multi <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/multiome_geo/multi_ref.rds")

# Load metadata
meta <- read.table("data/multiome_geo/metadata_csv/h5ad_file_metadata.csv", sep = ",", header = TRUE)
meta.subset <- meta %>% 
  group_by(cell_type) %>% 
  slice(1:2000)

table(meta.subset$cell_type)

# Select for some donor barcodes
obj.multi$BC <- rownames(obj.multi@meta.data)
obj.multi$to_keep <- obj.multi$BC %in% meta.subset$BC
table(obj.multi$to_keep)


# Subset
obj.multi.subset <- subset(obj.multi, subset = to_keep)


# Save
saveRDS(obj.multi.subset, "~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/multiome_geo/multi_ref_2000.rds")
