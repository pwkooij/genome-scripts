rm(list = ls())#clean workspace
setwd("/PATH/TO/WORK_DIRECTORY")
getwd()
dev.off()

library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)

## load funannotate annotations file
ann <- read.delim("/PATH/TO/FUNANNOTATE/annotate_results/genus_species_ID.annotations.txt", sep="\t", quote="", stringsAsFactors=FALSE)

## create COG subset ####
cog_df <- ann %>%
  filter(COG != "") %>%
  separate_rows(COG, sep="[;]") %>%
  count(COG, sort=F) %>%
  mutate(COG = str_replace(COG, "^[A-Z]:", ""))

# Plot
p_cog <- ggplot(cog_df, aes(x=factor(COG, levels = sort(unique(COG), decreasing = TRUE)), y=n, fill=COG)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=scales::hue_pal()(nrow(cog_df))) +
  theme_minimal(base_size=12) +
  labs(
    title="COG Functional Category Distribution",
    x="COG Functional Category",
    y="Number of Genes"
  ) +
  geom_text(aes(label=n), hjust=-0.1, size=3.5) +
  theme(
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    axis.text.y = element_text(size=10),
    legend.position = "none"
  )

p_cog

# Save
ggsave("COG_category_distribution.pdf", plot = p_cog, width=12, height=6, dpi=300)

## create CAZyme subset ####
cazyme_df <- ann %>%
  filter(CAZyme != "") %>%
  separate_rows(CAZyme, sep="[;]") %>%
  mutate(CAZyme = str_extract(CAZyme, "^[A-Z]+\\d+")) %>%
  filter(!is.na(CAZyme)) %>%
  count(CAZyme, sort=T)

# Plot
p_cazyme <- ggplot(cazyme_df[1:10,], aes(x=reorder(CAZyme, n), y=n)) +
  geom_bar(stat="identity", fill="#1a9850") +
  coord_flip() +
  theme_minimal(base_size=12) +
  labs(
    title="Top 10 represented CAZymes",
    x="CAZyme Category",
    y="Number of Genes"
  ) +
  geom_text(aes(label=n), hjust=-0.1, size=3.5) +
  theme(
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    axis.text.y = element_text(size=10),
    legend.position = "none"
  )

p_cazyme

# Save
ggsave("CAZyme_category_distribution.pdf", plot = p_cazyme, width=12, height=6, dpi=300)

# Create summarized CAZyme categories
cazyme_summary <- cazyme_df %>%
  mutate(
    Class = case_when(
      str_starts(CAZyme, "AA")  ~ "Auxiliary Activities",
      str_starts(CAZyme, "GH")  ~ "Glycoside Hydrolases",
      str_starts(CAZyme, "CBM") ~ "Carbohydrate-Binding Molecules",
      str_starts(CAZyme, "GT")  ~ "Glycosyl Transferases",
      str_starts(CAZyme, "CE")  ~ "Carbohydrate Esterases",
      str_starts(CAZyme, "PL")  ~ "Polysaccharide Lyases",
      TRUE                      ~ "Other"
    )
  ) %>%
  group_by(Class) %>%
  summarise(Total = sum(n)) %>%
  arrange(desc(Total))

# Plot
p_cazyme_summary <- ggplot(cazyme_summary, aes(x=reorder(Class, Total), y=Total, fill=Class)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values = scales::hue_pal()(nrow(cazyme_summary))) +
  theme_minimal(base_size=12) +
  labs(
    title="CAZyme Class Distribution",
    x="CAZyme Class",
    y="Number of Genes"
  ) +
  geom_text(aes(label=Total), hjust=-0.1, size=3.5) +
  theme(
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    axis.text.y = element_text(size=10),
    legend.position = "none"
  )

p_cazyme_summary

# Save
ggsave("CAZyme_class_summary.pdf", plot = p_cazyme_summary, width=12, height=6, dpi=300)

## Gene Ontologies GOs ####
# Extract and clean GO terms
go_df <- ann %>%
  filter(GO.Terms != "") %>%
  separate_rows(GO.Terms, sep="[;]") %>%
  mutate(
    Namespace = str_extract(GO.Terms, "GO_\\w+"),
    Namespace = str_remove(Namespace, "GO_"),
    Namespace = str_replace(Namespace, "component", "Cellular component"), # NEW
    Namespace = str_replace(Namespace, "function", "Molecular function"), # NEW
    Namespace = str_replace(Namespace, "process", "Biological process"), # NEW
    Term = str_extract(GO.Terms, "(?<=- ).*?(?= \\[)"),
    Term = str_trim(Term)
  ) %>%
  group_by(Namespace, Term) %>%
  summarise(Count = n(), .groups = "drop")

# Pick top 20 representative GO terms per namespace (adjust n = 10 or 15 as needed)
go_top <- go_df %>%
  group_by(Namespace) %>%
  slice_max(order_by = Count, n = 20) %>%
  arrange(Namespace, desc(Count))

# Assign colors for each namespace
namespace_colors <- c(
  "Molecular function" = "#1b9e77",   # green
  "Cellular component" = "#d95f02",  # orange
  "Biological process" = "#7570b3"     # blue
)

# Plot similar to the example figure
go_plot <- ggplot(go_top, aes(x = Count, y = reorder(Term, Count), fill = Namespace)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.2, color = "navy") +
  facet_wrap(~Namespace, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_fill_manual(values = namespace_colors, guide = "none") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12, hjust = 0),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Gene Ontology Summary",
    x = "Number of Genes"
  ) +
  coord_cartesian(xlim = c(0, max(go_top$Count) * 1.1))

go_plot

# Save
ggsave("GO_terms_summary.pdf", plot = go_plot, width=12, height=6, dpi=300)

## MEROPS peptidase ####
# Extract catalytic types
merops_df <- ann %>%
  filter(Protease != "") %>%
  separate_rows(Protease, sep="[;]") %>%
  mutate(
    Protease = str_trim(Protease),
    Type = str_sub(Protease, 1, 1)
  ) %>%
  mutate(
    Type = case_when(
      Type == "A" ~ "Aspartic",
      Type == "C" ~ "Cysteine",
      Type == "G" ~ "Glutamic",
      Type == "I" ~ "Inhibitor",
      Type == "M" ~ "Metallo",
      Type == "N" ~ "Asparagine",
      Type == "P" ~ "Mixed",
      Type == "S" ~ "Serine",
      Type == "T" ~ "Threonine",
      Type == "U" ~ "Unknown",
      Type == "X" ~ "Compound",
      TRUE ~ "Other"
    )
  ) %>%
  count(Type, sort = TRUE)

# Plot MEROPS summary
merops_plot <- ggplot(merops_df, aes(x = reorder(Type, n), y = n, fill = Type)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = scales::hue_pal()(nrow(merops_df))) +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "MEROPS Protease Distribution by Catalytic Type",
    x = "Catalytic Type",
    y = "Number of Genes"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  ) 

merops_plot

# Save
ggsave("MEROPS_summary.pdf", plot = merops_plot, width=12, height=6, dpi=300)

## PFAM clan summary ####
# parse PFAM entries
pfam_df <- ann %>%
  filter(PFAM != "") %>%
  separate_rows(PFAM, sep="[;]") %>%
  mutate(
    Pfam = str_trim(PFAM),
    Pfam_acc = str_extract(Pfam, "PF\\d+"),
    Pfam_acc = ifelse(is.na(Pfam_acc), Pfam, Pfam_acc)
    ) %>%
  filter(!is.na(Pfam_acc))

# load PFAM clan mapping table
pfam_clans <- read.delim("/mnt/databasedisk/Pfam-A.clans.tsv", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(pfam_clans) <- c("Pfam_acc", "Clan_acc", "Clan_name", "Pfam_name", "Pfam_description")

pfam_clans <- pfam_clans %>%
  filter(!is.na(Pfam_acc), Pfam_acc != "") %>%
  distinct(Pfam_acc, .keep_all = TRUE)

# join pfam_df with clan mapping table
pfam_with_clan <- pfam_df %>%
  left_join(pfam_clans, by = "Pfam_acc") %>%
  mutate(
    Clan_acc = ifelse(is.na(Clan_acc) | Clan_acc == "", "Unclanned", Clan_acc),
    Clan_name = ifelse(is.na(Clan_name) | Clan_name == "", "Unclanned", Clan_name),
    Clan_label = ifelse(Clan_acc == "Unclanned", "Unclanned", paste0(Clan_acc, ":", Clan_name))
  )

# create top 20 PFAM clan table
pfam_top20 <- pfam_with_clan %>%
  count(Clan_label, sort = TRUE) %>%
  slice_max(n, n = 20)

# Create a color palette similar to the MEROPS one
pfam_top20 <- pfam_top20 %>%
  mutate(fill_group = ifelse(Clan_label == "Unclanned", "Unclanned", Clan_label))

# Unique non-unclanned clans in plotting order
other_groups <- pfam_top20 %>%
  filter(fill_group != "Unclanned") %>%
  pull(fill_group) %>%
  unique()

# Generate pastel rainbow colors
n_colors <- length(other_groups)
rainbow_colors <- hue_pal()(n_colors)

# Name the colors so ggplot matches them correctly
named_colors <- setNames(rainbow_colors, other_groups)

# Add the grey color for "Unclanned"
color_palette <- c("Unclanned" = "#BEBEBE", named_colors)

# Plot
pfam_plot <- ggplot(pfam_top20, aes(x = reorder(Clan_label, n), y = n, fill = fill_group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = color_palette) +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top 20 PFAM Clans",
    x = "PFAM Clan",
    y = "Number of Genes"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

pfam_plot

# Save
ggsave("PFAM_summary.pdf", plot = pfam_plot, width=12, height=6, dpi=300)

## antiSMASH cluster summary ####
# Fast antiSMASH parser for region-level products
parse_antismash_products <- function(gbk_file) {
  products <- c()
  current_section <- NULL
  
  con <- file(gbk_file, open = "r")
  while (TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break  # End of file
    
    # Detect section (region / protocluster / proto_core etc.)
    if (grepl("^\\s{5}\\w", line)) {
      current_section <- strsplit(trimws(line), "\\s+")[[1]][1]
    }
    
    # Capture product only if inside 'region'
    if (!is.null(current_section) && current_section == "region" && grepl("/product=", line)) {
      product <- sub('.*\\/product="([^"]+)".*', "\\1", line)
      products <- c(products, product)
    }
  }
  close(con)
  
  products_df <- as.data.frame(table(products))
  colnames(products_df) <- c("Product", "Count")
  products_df <- products_df[order(products_df$Count, decreasing = TRUE), ]
  return(products_df)
}

# Example use
products_df <- parse_antismash_products("/PATH/TO/antismash/antismash_ID.gbk")
head(products_df)

antismash_plot <- ggplot(products_df, aes(x = reorder(Product, Count), y = Count, fill = Product)) +
  geom_col(color = "black", linewidth = 0.3) +
  coord_flip() +
  geom_text(aes(label = Count), hjust = -0.2, size = 3) +
  scale_fill_manual(values = scales::hue_pal()(nrow(products_df))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "antiSMASH Secondary Metabolite Clusters by Product Type",
    x = "Cluster Type",
    y = "Number of Regions"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

antismash_plot

# Save
ggsave("antiSMASH_summary.pdf", plot = antismash_plot, width=12, height=6, dpi=300)
