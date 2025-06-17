library(tidyr)
library(tidyverse)
library(vegan)
library(cluster)
library(fastDummies)
library(tibble)
library(ggplot2)
library(SYNCSA)
library(patchwork)
library(dplyr)

# Bring in data

birds <- read.csv("Master_DB_Birds.csv")
name <- read.csv("tran.birdlist.alldetails.csv")

## Select summer seasons
birds <- birds %>%
  filter(period %in% c("1", "3", "4"))

# merge names
birds <- birds %>%
  rename(common_name = species)

birds <- birds %>%
  left_join(name %>% select(common_name, species, genus, family, order), by = "common_name")

## concat column and add _ to name
birds$tran_id <- paste(birds$transect, birds$water, sep = "_")
birds$species <- gsub(" ", "_", birds$species)  # Replace spaces with underscores


## Create p/a matrix (1-0)
birds$presence <- ifelse(birds$count > 0, 1, 0)
presence_absence_matrix <- xtabs(presence ~ tran_id + species, data = birds)

# Create the abundance matrix using actual counts
abundance_matrix <- xtabs(count ~ tran_id + species, data = birds)

# Convert the matrix to a data frame if needed for further processing or output
abundance_df <- as.data.frame.matrix(abundance_matrix)


# Convert counts to 1 (presence) and 0 (absence)
presence_absence_matrix[presence_absence_matrix > 0] <- 1

# View the matrix
print(presence_absence_matrix)
write.csv(as.data.frame.matrix(presence_absence_matrix), "presence_absence_matrix.csv", row.names = TRUE)
write.csv(as.data.frame.matrix(abundance_matrix), "abundance_matrix.csv", row.names = TRUE)
pa <- read.csv("presence_absence_matrix.csv", row.names = 1, check.names = FALSE)
ab <- read.csv("abundance_matrix.csv", row.names = 1, check.names = FALSE)


## AvTD - average taxonomic distance (Distinctiveness)

nodes_unique <- birds %>%
  select(species, genus, family, order) %>%
  distinct()
nodes_unique <- na.omit(nodes_unique)

avtd_matrix <- nodes_unique %>%
  filter(species %in% colnames(pa)) %>%
  arrange(factor(species, levels = colnames(pa)))

# Set species as row names
rownames(avtd_matrix) <- avtd_matrix$species
avtd_matrix <- avtd_matrix[, -1]
avtd_dist_matrix <- taxa2dist(as.matrix(avtd_matrix))
avtd_dist_matrix <- as.matrix(avtd_dist_matrix)

# Check if matrix is symmetric
all.equal(avtd_dist_matrix, t(avtd_dist_matrix))

# Check if diagonal is all zeros
all(diag(avtd_dist_matrix) == 0)


# Calculate Average Taxonomic Distinctiveness (ATD) for each sample
avtd_results <- apply(pa, 1, function(sample) {
  present_species <- which(sample == 1)
  
  if (length(present_species) > 1) {
    # Extract the distances for the present species
    distances <- avtd_dist_matrix[present_species, present_species]
    
    # Calculate the average taxonomic distance
    avtd <- mean(distances[upper.tri(distances)], na.rm = TRUE)  # Use upper.tri to avoid redundancy
    return(avtd)
  } else {
    return(NA)  # Return NA if less than 2 species are present
  }
})


# dataframe
avtd_df <- data.frame(tran_id = rownames(pa), AvTD = avtd_results)
# Merge with metadata

# Pull transect conditions
unique_transects <- birds %>%
  distinct(tran_id, .keep_all = TRUE) %>%
  select(tran_id, water)
avtd_df <- left_join(avtd_df, unique_transects, by = "tran_id")

p1 <- ggplot(avtd_df, aes(x = water, y = AvTD, fill = water)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  
  geom_point(alpha = 0.5, size = 2.5, position = position_jitter(width = 0.03)) +
  scale_fill_manual(values = c("Yes" = "lightblue", "No" = "lightgreen")) +
  labs(x = NULL, y = "AvTD") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 14,margin = margin(r = 10)),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0),
    legend.position = "none"
  )

print(p1)


ggsave('summer_avTD_boxplot.png', plot = p1, height = 4, width = 4, dpi = 300, bg = "white")

## Final format

# Ensure the factor levels are ordered with "Yes" first
avtd_df$water <- factor(avtd_df$water, levels = c("Yes", "No"), labels = c("With Water", "Without Water"))

# Create pair_id by taking the integer part of tran_id and grouping each pair
avtd_df <- avtd_df %>%
  mutate(pair_id = ceiling(as.numeric(gsub("_.*", "", tran_id)) / 2))

# Now plot with connecting lines
p1_s <- ggplot(avtd_df, aes(x = water, y = AvTD, fill = water)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_line(aes(group = pair_id), color = "gray40", linewidth = 0.5, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.3) +
  scale_fill_manual(values = c("With Water" = "#2A758E", "Without Water" = "#65A73D")) +
  labs(x = NULL, y = "AvTD") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0),
    legend.position = "none"
  )
p1_s
ggsave("summer_avtd.png", plot = p1_s, width = 6, height = 8, units = "in", dpi = 300)

## Trait dissimilarity matrix for FD and Rao's quadratic entropy ----

# Step 2: Select relevant functional trait columns
trait_data <- name %>%
  select(species, Mass, Trophic.Niche, Primary.Lifestyle, size_class) %>% 
  mutate(Mass = log(Mass + 1))
trait_data$species <- gsub(" ", "_", trait_data$species)


# Convert relevant columns to factors
# trait_data$Mass_bins <- factor(trait_matrix$Mass_bins, ordered = TRUE)
trait_data$Trophic.Niche <- as.factor(trait_data$Trophic.Niche)
trait_data$Primary.Lifestyle <- as.factor(trait_data$Primary.Lifestyle)
trait_data$size_class <- as.factor(trait_data$size_class)

# Check mass is numeric
trait_data$Mass <- as.numeric(trait_data$Mass)

# Calculate the Gower similarity matrix
gower_dissimilarity <- daisy(trait_data %>% select(-species, -Mass), metric = "gower")
gower_dissimilarity_matrix <- as.matrix(gower_dissimilarity)

rownames(gower_dissimilarity_matrix) <- trait_data$species
colnames(gower_dissimilarity_matrix) <- trait_data$species

## AvFD Average Functional Distance ----

avfd_results <- apply(pa, 1, function(sample) {
  present_species <- which(sample == 1)
  present_species_names <- colnames(pa)[present_species]  # Get species names
  
  # Check if at least 2 species are present
  if (length(present_species_names) > 1) {
    # Extract similarity values for present species
    similarities <- gower_dissimilarity_matrix[present_species_names, present_species_names]
    avfd <- mean(similarities[upper.tri(similarities)], na.rm = TRUE)  # Use upper triangle to avoid redundancy
    return(avfd)
  } else {
    return(NA)  # Return NA if fewer than 2 species
  }
})

# arrange df
avfd_df <- data.frame(sample = rownames(pa), AvFD = avfd_results)
avfd_df <- avfd_df %>%
  rename(tran_id = sample)
avfd_df <- left_join(avfd_df, unique_transects, by = "tran_id")

# plot
p2 <- ggplot(avfd_df, aes(x = water, y = AvFD, fill = water)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  
  geom_point(alpha = 0.5, size = 2.5, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Yes" = "lightblue", "No" = "lightgreen")) +
  labs(x = NULL, y = "AvFD") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0),
    legend.position = "none"
  )

print(p2)

ggsave('average_functional_diversity.png', plot = avFD_boxplot, height = 4, width = 4, dpi = 300, bg = "white")


#### rao with SYNCSA
## Rao's quadratic entropy ----
# Ensure species are row names in traits
trait_data <- name %>%
  select(species, Mass, Trophic.Niche, Primary.Lifestyle, size_class) %>% 
  mutate(Mass = log(Mass + 1))
trait_data$species <- gsub(" ", "_", trait_data$species)

traits <- trait_data %>%
  column_to_rownames(var = "species")

traits <- traits %>%
  mutate(across(where(is.character), as.factor))

# Step 3: Calculate RaoQ for filtered samples
rao_results <- rao.diversity(ab, traits, ord = "metric", standardize = TRUE)
rao_values <- rao_results$FunRao  # Extract RaoQ values

# Step 4: Convert rao_values to a data frame
rao_df <- data.frame(sample = names(rao_values), RaoQ = rao_values, row.names = NULL)

# Step 5: Remove zero values
rao_df_filtered <- rao_df %>%
  filter(RaoQ > 0)

rao_df_filtered <- rao_df_filtered %>%
  rename(tran_id = sample)

# Step 6: Join with metadata
rao_df_filtered <- left_join(rao_df_filtered, unique_transects, by = "tran_id")

# Step 8: Create the boxplot
p3 <- ggplot(rao_df_filtered, aes(x = water, y = RaoQ, fill = water)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  
  geom_point(alpha = 0.5, size = 2.5, position = position_jitter(width = 0.03)) +
  scale_fill_manual(values = c("Yes" = "lightblue", "No" = "lightgreen")) +
  labs(x = NULL, y = "Functional Diversity (RaoQ)") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0),
    legend.position = "none"
  )

print(p3)

ggsave('rao_functional_diversity.png', plot = p3, height = 4, width = 4, dpi = 300, bg = "white")

## Redundancy # Extract functional redundancy values
redundancy_values <- rao_results$FunRedundancy

# Turn into a dataframe
redundancy_df <- data.frame(
  tran_id = names(redundancy_values),
  FunctionalRedundancy = redundancy_values
)

# Join with metadata (if needed)
redundancy_df <- left_join(redundancy_df, unique_transects, by = "tran_id")

# Plot
p4 <- ggplot(redundancy_df, aes(x = water, y = FunctionalRedundancy, fill = water)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  
  geom_point(alpha = 0.5, size = 2.5, position = position_jitter(width = 0.03)) +
  scale_fill_manual(values = c("Yes" = "lightblue", "No" = "lightgreen")) +
  labs(x = NULL, y = "Functional Redundancy") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0)
  )

print(p4)

ggsave('redundancy_functional_diversity.png', plot = p4, height = 4, width = 4, dpi = 300, bg = "white")

## put the plots in a row

library(scales)  #
library(patchwork)
# Example for each plot
p1 <- p1 + scale_y_continuous(labels = number_format(accuracy = 0.01))
p3 <- p3 + scale_y_continuous(labels = number_format(accuracy = 0.01))
p4 <- p4 + scale_y_continuous(labels = number_format(accuracy = 0.01))


combined_plot <- p1 + p3 + p4 
plot_layout(ncol = 3)
combined_plot
ggsave('combined_summer.png', plot = combined_plot, height = 4, width = 16, dpi = 300, bg = "white")

# __________________________________________________________________

# Diagnostic plots
avtd_df$AvTD <- as.numeric(avtd_df$AvTD)
hist(avtd_df$AvTD[avtd_df$water == "Yes"])
hist(avtd_df$AvTD[avtd_df$water == "No"])
boxplot(avtd_df$AvTD, water = "Yes")
boxplot(log(avtd_df$AvTD) ~ redundancy_df$water)
leveneTest(avtd_df$AvTD ~ water, data = redundancy_df)

hist(rao_df_filtered$RaoQ[rao_df_filtered$water == "Yes"])
hist(rao_df_filtered$RaoQ[rao_df_filtered$water == "No"])
boxplot(rao_df_filtered$RaoQ ~ rao_df_filtered$water)
boxplot(log(rao_df_filtered$RaoQ) ~ rao_df_filtered$water)

hist(redundancy_df$FunctionalRedundancy[redundancy_df$water == "Yes"])
hist(redundancy_df$FunctionalRedundancy[redundancy_df$water == "No"])
boxplot(redundancy_df$FunctionalRedundancy ~ redundancy_df$water)
boxplot(log(redundancy_df$FunctionalRedundancy) ~ redundancy_df$water)

# ______________________________________________________________________________
# Paired test - Taxanomic distinctiveness
avtd_df <- avtd_df %>%
  mutate(pair_id = as.numeric(sub("_.*", "", tran_id)),  # Extract numeric part of ID
         pair_id = ceiling(pair_id / 2),  # Create correct pair numbering (1 for 1_Yes and 2_No, 2 for 3_Yes and 4_No, etc.)
         water_status = ifelse(grepl("_Yes", tran_id), "Yes", "No"))  # Identify Yes/No

# Step 2: Pivot to align Yes and No in the same row
paired_data <- avtd_df %>%
  select(pair_id, water_status, AvTD) %>%  # Keep relevant columns
  pivot_wider(names_from = water_status, values_from = AvTD, names_prefix = "AvTD_") %>%
  arrange(pair_id)  # Ensure correct ordering

wilcox.test(paired_data$AvTD_Yes, paired_data$AvTD_No, paired = TRUE)

# Paired test - Rao functional diversity
rao_df_filtered <- rao_df_filtered %>%
  mutate(pair_id = as.numeric(sub("_.*", "", tran_id)),  # Extract numeric part of ID
         pair_id = ceiling(pair_id / 2),  # Create correct pair numbering (1 for 1_Yes and 2_No, 2 for 3_Yes and 4_No, etc.)
         water_status = ifelse(grepl("_Yes", tran_id), "Yes", "No"))  # Identify Yes/No

# Paired test - Taxanomic distinctiveness
# Step 2: Pivot to align Yes and No in the same row
paired_data_rao <- rao_df_filtered %>%
  select(pair_id, water_status, RaoQ) %>%  # Keep relevant columns
  pivot_wider(names_from = water_status, values_from = RaoQ, names_prefix = "RaoQ_") %>%
  arrange(pair_id)  # Ensure correct ordering

wilcox.test(paired_data_rao$RaoQ_Yes, paired_data_rao$RaoQ_No, paired = TRUE)

# Paired test - Functional redundancy
redundancy_df <- redundancy_df %>%
  mutate(pair_id = as.numeric(sub("_.*", "", tran_id)),  # Extract numeric part of ID
         pair_id = ceiling(pair_id / 2),  # Create correct pair numbering (1 for 1_Yes and 2_No, 2 for 3_Yes and 4_No, etc.)
         water_status = ifelse(grepl("_Yes", tran_id), "Yes", "No"))  # Identify Yes/No

# Paired test - Taxanomic distinctiveness
# Step 2: Pivot to align Yes and No in the same row
paired_data_red <- redundancy_df %>%
  select(pair_id, water_status, FunctionalRedundancy) %>%  # Keep relevant columns
  pivot_wider(names_from = water_status, values_from = FunctionalRedundancy, names_prefix = "Fun_") %>%
  arrange(pair_id)  # Ensure correct ordering

wilcox.test(paired_data_red$Fun_Yes, paired_data_red$Fun_No, paired = TRUE)

# ______________________________________________________________________________
