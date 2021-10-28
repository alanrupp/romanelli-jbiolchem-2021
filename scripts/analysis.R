# - Analysis of RNA-seq data from BAT after UCP1 KO ---------------------------
library(dplyr)
library(ggplot2)
library(stringr)
if (!dir.exists("results")) dir.create("results")

# - Read in data --------------------------------------------------------------
metadata <- read.csv("data/samples.csv")
genes <- read.csv("data/genes.csv") %>%
  filter(!duplicated(gene_name)) %>%
  select(starts_with("gene"))

# read in count data
read_counts <- function(sample, column = "V4") {
  file <- paste0("alignments/", sample, "/ReadsPerGene.out.tab")
  df <- read.table(file, skip = 4, sep = "\t")
  df <- select(df, V1, !!sym(column))
  colnames(df) <- c("gene_id", sample)
  df
}
counts <- lapply(metadata$Sample, read_counts) %>% bind_cols()

counts <- select(counts, "gene_id...1", starts_with("1"))
counts <- rename(counts, "gene_id" = `gene_id...1`)

# - Quality control -----------------------------------------------------------
colors <- c("white", "#531b93")

# library size
plot_library_size <- function(counts) {
  size <- colSums(select(counts, -gene_id))
  size <- as.data.frame(size) %>%
    tibble::rownames_to_column("Sample") %>%
    mutate(size = size / 10^6) %>%
    left_join(metadata, by = "Sample") %>%
    mutate(Sample = factor(Sample, levels = metadata$Sample))
  # plot
  ggplot(size, aes(x = Sample, y = size, fill = Treatment)) +
    geom_hline(aes(yintercept = 30), linetype = "dashed") +
    geom_col(color = "black") +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = round(size, 0), y = size + 3)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, max(size$size) + 5)) +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("Library size (M)")
}
p <- plot_library_size(counts)
ggsave("results/library_size.png", p, width = 4, height = 2, units = "in", dpi = 600)

# comparing similarity of WT and Ucp1 KO samples (Ucp1 should be down)
# make counts into cpm
make_cpm <- function(counts, log2 = FALSE) {
  if (log2) {
    df <- mutate_if(counts, is.numeric, ~ log2(.x / sum(.x) * 10^6 + 1))
  } else {
    df <- mutate_if(counts, is.numeric, ~ .x / sum(.x) * 10^6)
  }
  df
}

#' plot a gene's CPM s a violin for each treatment
plot_gene <- function(gene) {
  df <- make_cpm(counts) %>%
    left_join(genes, by = "gene_id") %>%
    filter(gene_name == gene) %>%
    tidyr::pivot_longer(-starts_with("gene"), names_to = "Sample", values_to = "cpm") %>%
    left_join(metadata, by = "Sample")
  # plot
  xlabels <- c("BAd-CRISPR\nControl", "BAd-CRISPR\nUcp1")
  ggplot(df, aes(x = Treatment, y = cpm)) +
    geom_violin(aes(fill = Treatment), show.legend = FALSE, scale = "width") +
    scale_fill_manual(values = colors, guide = "none") +
    geom_jitter(show.legend = FALSE) +
    theme_classic() +
    scale_x_discrete(labels = xlabels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    ylab("CPM") +
    ggtitle(gene) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 10, face = "italic"),
          axis.text.x = element_text(size = 9, color = "black", face = "bold"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          axis.title = element_text(size = 9, face = "bold"),
          axis.text.y = element_text(size = 8, color = "black", face = "bold")) +
    xlab(NULL)
}

p <- plot_gene("Ucp1")
ggsave("results/Ucp1_CPM.png", p, width = 3.3, height = 2.5, units = "in", dpi = 600)

# sample similarity (PCA embedding)
find_expressed <- function(counts, threshold = 1, samples = 3) {
  cpm <- make_cpm(counts, log2 = FALSE)
  above_threshold <- mutate_if(cpm, is.numeric, ~ .x >= threshold)
  total_samples <- rowSums(select(above_threshold, -gene_id))
  return(cpm$gene_id[total_samples >= samples])
}
expressed_genes <- find_expressed(counts)
cpm <- make_cpm(counts, log2 = TRUE) %>% filter(gene_id %in% expressed_genes)
pca <- prcomp(t(cpm[, 2:ncol(cpm)]))

df <- pca$x[, c("PC1", "PC2")] %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(metadata, by = "Sample")

# get % variance explained
variance <- pca$sdev^2
percent_var <- (variance / sum(variance) * 100) %>% round(0)
names(percent_var) <- colnames(pca$x)

p <- ggplot(df, aes(x = PC1, y = PC2, fill = Treatment)) +
  geom_point(shape = 21) +
  ggrepel::geom_text_repel(aes(label = Sample), show.legend = FALSE) +
  theme_classic() +
  xlab(paste0("PC1 (", percent_var["PC1"], "% of variance)")) +
  ylab(paste0("PC2 (", percent_var["PC2"], "% of variance)")) +
  scale_fill_manual(values = colors)
ggsave("results/PCA.png", p, width = 4, height = 3, units = "in", dpi = 600)

# - DE analysis ---------------------------------------------------------------
get_de <- function() {
  counts <- as.data.frame(counts) %>% tibble::column_to_rownames("gene_id")
  samples <- as.data.frame(metadata) %>% tibble::column_to_rownames("Sample")
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, samples, ~ Treatment) %>%
    DESeq2::DESeq()
  de <- DESeq2::results(dds) %>% as.data.frame()
  cpm <- DESeq2::fpm(dds)
  list("dds" = dds, "de" = de, "cpm" = cpm)
}
de <- get_de()

plot_zscores <- function(gene_names) {
  df <- de$cpm %>% as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(genes, by = "gene_id") %>%
    filter(gene_name %in% gene_names) %>%
    tidyr::pivot_longer(-starts_with("gene"), names_to = "Sample", values_to = "CPM") %>%
    mutate(gene_name = factor(gene_name, levels = gene_names[length(gene_names):1])) %>%
    left_join(metadata, by = "Sample") %>%
    arrange(Treatment) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))
  # get z-scores
  df <- df %>% group_by(gene_name) %>% mutate("z" = (CPM-mean(CPM))/sd(CPM))
  # plot
  ggplot(df, aes(x = Sample, y = gene_name, fill = z)) +
    geom_tile() +
    scale_fill_gradient(low = colors[1], high = colors[2]) +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(hjust = 1, size = 6, face = "italic"),
          legend.key.width = unit(0.02, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
}

top_genes <- list(
  "Perixisome" = c("Pex16", "Nudt7", "Abcd2", "Acot8", "Abcd3", "Pex19",
                   "Acsl5", "Hsd17b4", "Pex13", "Pex14", "Acot4"),
  "Lipid Metabolism" = c("Plin2", "Fabp4", "Gpd1", "Agpat2", "Agpat3", "Acot11",
                         "Me1", "Acot2", "Pla2g2e", "Hsdl2", "Hmgcs1", "Elovl5",
                         "Hccs", "Pdhb", "Hsd17b7", "Elovl3", "Scd3"),
  "Mitochondria" = c("mt-Nd4", "mt-Nd1", "mt-Rnr1", "mt-Co1", "mt-Rnr2",
                     "mt-Nd3", "mt-Atp6", "mt-Nd6", "Gcdh"),
  "Translation" = c("Cars", "Nars", "Yars", "Wars", "Rps6ka1", "Tars", "Sars",
                    "Gars", "Iars", "Eprs"),
  "Proteasome" = c("Psmd12", "Psmd11", "Psmd14", "Psmd13", "Psma5", "Psmb7",
                   "Psma6", "Psmd6", "Psmb5", "Psmd7")
)

p <- plot_zscores(unlist(top_genes))
ggsave("results/de_zscores.png", width = 3.6, height = 4.1, units = "in", dpi = 600)


# off-target genes
off_target <- c("Ucp1", "Ypel2", "Katnal1", "Ago3", "Adamts14", "Neu2", "Rftn1")
p <- plot_zscores(off_target)
ggsave("results/offtarget_zscores.png", width = 2, height = 2, units = "in", dpi = 600)

# - GSEA ----------------------------------------------------------------------
# order genes by DESeq2 stat
ranks <- de$de$stat
names(ranks) <- rownames(de$de)
ranks <- ranks[!is.na(de$de$padj)]
ranks <- sort(ranks)

# GSEAPreranked data
rank_df <- data.frame("gene" = names(ranks), "rank" = ranks)
write.table(rank_df, "results/ranks.rnk", sep = "\t", row.names = FALSE, quote = FALSE)

# make CLS file to describe the classes
line1 <- c(nrow(metadata), length(unique(metadata$Treatment)), 1) %>%
  paste(collapse = " ")
line2 <- c("#", str_remove(sort(unique(metadata$Treatment)), " ")) %>%
  paste(collapse = " ")
line3 <- mutate(metadata, Treatment = ifelse(Treatment == "Control", 0, 1)) %>%
  pull(Treatment) %>%
  paste(collapse = " ")
file_connection <- file("results/groups.cls")
writeLines(c(line1, line2, line3), file_connection)
close(file_connection)

# !!! GSEA is run offline with the Broad software

# - Save output for GEO -------------------------------------------------------
rename(counts, "id" = gene_id) %>% write.csv("results/STAR_counts.csv", row.names = FALSE)
