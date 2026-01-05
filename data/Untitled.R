library(Seurat)
library(Matrix)
getwd()
setwd("/Users/zz005/Documents/github/docker/my_docker_3.0/data/")

data_dir <- "/Users/zz005/Documents/github/docker/my_docker_3.0/data/"

expression_matrix <- Read10X(data.dir = data_dir)

# 3. Create a Seurat Object
# We filter slightly here to ensure the object is valid (at least 3 cells per gene)
seurat_obj <- CreateSeuratObject(
  counts = expression_matrix, 
  project = "test_run"
)

# 4. Save the object as an RDS file
output_file <- paste0(data_dir, "matrix.rds")
print(paste("Saving RDS to:", output_file))
saveRDS(seurat_obj, file = output_file)

print("Done!")
