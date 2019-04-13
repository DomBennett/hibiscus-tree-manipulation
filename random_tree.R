# Libs ----
library(ape) # TODO: switch to treeman

# Functions ----
spnms_get <- function(tree) {
  tip_labels <- tree$tip.label
  res <- sub(pattern = '_copy_?[0-9]+', replacement = '',
             x = tip_labels)
  res <- sub(pattern = "_MAP(P_)?[0-9]+", replacement = '', x = res)
  res <- sub(pattern = "(\\(reversed\\))", replacement = '', x = res)
  res <- sub(pattern = "_blast.*$", replacement = '', x = res)
  res <- sub(pattern = "_copy$", replacement = '', x = res)
  res
}
random_species_tree_get <- function(tree, n_tips_per_species) {
  spp <- spnms_get(tree)
  unique_spp <- unique(spp)
  tips_to_keep <- NULL
  for (sp in unique_spp) {
    pssbls <- tree$tip.label[spp == sp]
    nmax <- ifelse(length(pssbls) >= n_tips_per_species, n_tips_per_species,
                   length(pssbls))
    tips_to_keep <- c(tips_to_keep, sample(x = pssbls, size = nmax))
  }
  tips_to_drop <- tree$tip.label[!tree$tip.label %in% tips_to_keep]
  res <- ape::drop.tip(phy = tree, tip = tips_to_drop)
  res$tip.label <- spnms_get(res)
  res
}

# Vars ----
nrands <- 100  # Number of rand. trees per gene tree
output_file <- "random_trees.newick"
n_tips_per_species <- 2

# Data ----
gene_tree_files <- file.path('data',
                             list.files(path = 'data', pattern = '.newick'))

# Loop ----
if (file.exists(output_file)) file.remove(output_file)
counter <- 0
for (fl in gene_tree_files) {
  gene_tree <- ape::read.tree(file = fl)
  for (i in seq_len(nrands)) {
    rand_tree <- random_species_tree_get(tree = gene_tree, n_tips_per_species =
                                           n_tips_per_species)
    write.tree(phy = rand_tree, file = output_file, append = TRUE)
    counter <- counter + 1
  }
}

# all_trees <- list()
# for (fl in gene_tree_files) {
#   all_trees[[fl]] <- ape::read.tree(file = fl)
# }
# spp_nms <- lapply(X = all_trees, spnms_get)
# unique(unlist(spp_nms))
