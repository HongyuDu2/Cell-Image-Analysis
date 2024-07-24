# library
library(dplyr)
library(tidyr)
library(readxl)

# Import dataset
cell_segmentation_table <- read.csv("C:/Users/duho/Desktop/cell_segmentation_table.csv")
clustering_Table <- read_excel("C:/Users/duho/Desktop/Final_Clustering_Table.xlsx")
unidentified_cell <- read.csv("C:/Users/duho/Desktop/unidentified_cells_table.csv")
unidentified_clustering <- read_excel("C:/Users/duho/Desktop/Unidentified_Clustering_Table-2.xlsx")

cell_type <- c()
for (i in 1:nrow(cell_segmentation_table)){
  for (j in 1:200){
    if (cell_segmentation_table$cluster[i] == clustering_Table$cluster[j]){
      cell_type[i] = clustering_Table$Cell_Type_renew1[j]
      break
    }
  }
}
cell_segmentation_table$cell_type_renew1 <- cell_type

cell_type <- c()
for (i in 1:nrow(unidentified_cell)){
  for (j in 1:nrow(unidentified_clustering)){
    if (unidentified_cell$cluster_unidentified[i] == unidentified_clustering$cluster[j]){
      cell_type[i] = unidentified_clustering$revised_cell_type_1[j]
      break
    }
  }
}
unidentified_cell$cell_type_renew1 <- cell_type

cell_segmentation_table$cell_type_renew2 <- cell_segmentation_table$cell_type_renew1
cell_segmentation_table_x <- cell_segmentation_table %>% select('label', 'cell_type_renew2', 'fov')
unidentified_cell_x <- unidentified_cell %>% select('label' , 'cell_type_renew1', 'fov')
merged_table <- left_join(cell_segmentation_table_x, unidentified_cell_x, by = c("label", "fov"))
merged_table <- merged_table %>%
  mutate(cell_type_renew2 = if_else(!is.na(cell_type_renew1), cell_type_renew1, cell_type_renew2)) %>%
  select(label, celltype = cell_type_renew2)
cell_segmentation_table$cell_type_renew2 <- merged_table$celltype

write.csv(cell_segmentation_table, "C:/Users/duho/Desktop/cell_segmentation_table.csv", row.names = FALSE)
write.csv(unidentified_cell, "C:/Users/duho/Desktop/unidentified_cells_table.csv", row.names = FALSE)


# Calculate the cell_type amount in each fov
cell_counts <- cell_segmentation_table %>%
  group_by(fov, cell_type_renew2) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total amount of cells in each fov
total_counts <- cell_segmentation_table %>%
  group_by(fov) %>%
  summarise(total = n(), .groups = 'drop')

# Combine two dataset and calculate the percentages
cell_percentages <- cell_counts %>%
  left_join(total_counts, by = "fov") %>%
  mutate(percentage = (count / total) * 100) %>%
  select(fov, cell_type_renew2, percentage)

# Change to wide-form
cell_percentage_wide <- cell_percentages %>%
  pivot_wider(names_from = cell_type_renew2, values_from = percentage, values_fill = list(percentage = 0))

write.csv(cell_percentage_wide, "C:/Users/duho/Desktop/percentage_with_unidentified.csv", row.names = FALSE)


#percentage without unidentified

cell_segmentation_table_short <- cell_segmentation_table %>% filter(cell_type_renew2 != "Unidentified")

# Calculate the cell_type amount in each fov
cell_counts <- cell_segmentation_table_short %>%
  group_by(fov, cell_type_renew2) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total amount of cells in each fov
total_counts <- cell_segmentation_table_short %>%
  group_by(fov) %>%
  summarise(total = n(), .groups = 'drop')

# Combine two dataset and calculate the percentages
cell_percentages <- cell_counts %>%
  left_join(total_counts, by = "fov") %>%
  mutate(percentage = (count / total) * 100) %>%
  select(fov, cell_type_renew2, percentage)

# Change to wide-form
cell_percentage_wide <- cell_percentages %>%
  pivot_wider(names_from = cell_type_renew2, values_from = percentage, values_fill = list(percentage = 0))

write.csv(cell_percentage_wide, "C:/Users/duho/Desktop/percentage_without_unidentified.csv", row.names = FALSE)



# create the table
summary_table <- cell_segmentation_table %>%
  group_by(fovs, cell_type) %>%
  summarise(count = n()) %>%
  spread(key = cell_type, value = count, fill = 0)
# save the summary table
write.csv(summary_table, "C:/Users/duho/Desktop/fov_celltype_table.csv", row.names = FALSE)


# Created Unidentified Cells Table
unidentified_cells_table <- cell_segmentation_table[cell_segmentation_table$cell_type_renew1 == "Unidentified", ]
print(unidentified_cells)
