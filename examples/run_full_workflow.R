# =========================================================
# MicroSSNet: End-to-end example workflow
# This script demonstrates how to run:
# 1) SSN analysis
# 2) Aggregated network analysis
# 3) Bipartite network analysis
# from example input files to final outputs
# =========================================================

# -------------------------
# 0. Install dependencies
# -------------------------
install.packages("BiocManager")
BiocManager::install("limma")
install.packages("remotes")

remotes::install_github("zdk123/SpiecEasi")
remotes::install_github("TangZecheng622/MicroSSNet")

library(MicroSSNet)
library(phyloseq)

# -------------------------
# 1. Load example data
# -------------------------
abund <- read.csv("example_input_data.csv", row.names = 1, check.names = FALSE)
meta  <- read.csv("example_input_meta.csv", stringsAsFactors = FALSE)
tax   <- read.csv("example_tax_input_data.csv", row.names = 1, check.names = FALSE)
abund_fungi <- read.csv("example_fungi_input_data.csv", row.names = 1, check.names = FALSE)

# -------------------------
# 2. Build phyloseq object
# -------------------------
otu <- otu_table(as.matrix(abund), taxa_are_rows = TRUE)
sam <- sample_data(meta)
rownames(sam) <- sam$sample
phyloseq_test <- phyloseq(otu, sam)

# -------------------------
# 3. SSN analysis
# -------------------------
# option 1: abundance table + metadata
# group_df is required when table1 is a plain abundance table
ssn_test <- ssn_pipeline(
  table1 = abund,
  group_df = meta,
  ssn_method = "ssPCC",   # or "LIONESS-S", "LIONESS-D"
  vscol1 = "year",        # metadata column used for grouping
  control = "2009",       # baseline group
  vscol2 = "Warm"         # second metadata column for subgroup comparison
)

# option 2: phyloseq object
# if sample metadata are already stored in the phyloseq object,
# group_df can be set to NULL
ssn_test <- ssn_pipeline(
  table1 = phyloseq_test,
  group_df = NULL,
  ssn_method = "ssPCC",
  vscol1 = "year",
  control = "2009",
  vscol2 = "Warm"
)

# -------------------------
# 4. Aggregated network analysis
# Supports two input modes:
#   option 1: abundance table + metadata
#   option 2: phyloseq object
#
# tax is optional:
#   - if provided, node colors can be assigned according to a taxonomic column
#   - if omitted, nodes will be colored by network modules instead
# -------------------------

# option 1: phyloseq object
# group_df can be NULL if sample metadata are already stored in the phyloseq object
agg_net_test <- aggregation_netpipeline(
  table1 = phyloseq_test,
  group_df = NULL,
#  tax = tax,              # optional
#  tax_col = "Phylum",     # optional; used only when tax is provided
  vscol1 = "Warm",
  method = "spearman",
  r.threshold = 0.5,
  p.threshold = 0.5,
  step = 50,
  output_dir = "./example_output/aggregated_network_phyloseq"
)

# option 2: abundance table + metadata
agg_net_test <- aggregation_netpipeline(
  table1 = abund,
  group_df = meta,
  tax = tax,              # optional
  tax_col = "Phylum",     # optional; used only when tax is provided
  vscol1 = "Warm",
  method = "spearman",
  r.threshold = 0.5,
  p.threshold = 0.5,
  step = 50,
  output_dir = "./example_output/aggregated_network_table"
)

# -------------------------
# 5. Bipartite network analysis
# -------------------------
b_net <- bipartite_netpipeline(
  table1 = abund,
  table2 = abund_fungi,
  table1name = "bac",
  table2name = "fungi",
  group_df = meta,
  vscol1 = "Warm",
  method = "spearman",
  r.threshold = 0.5,
  p.threshold = 0.5,
  output_dir = "./example_output/bipartite_network"
)
