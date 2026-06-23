# ============================================================
# Colocalization analysis: DLB GWAS x cell-type-specific eQTLs
# Reproduces the coloc results reported in Table 1
# ============================================================

library(data.table)
library(coloc)

# ---- 1. Load GWAS summary statistics ----
gwas_dt <- fread("data/LBD_GWAS_SMR_formatted.txt")

# ---- 2. Load reusable analysis functions ----
source("scripts/coloc_functions.R")

# ---- 3. Discover all gene x cell-type combinations available ----
esd_dir <- "data/eqtl_esd"
esd_files <- list.files(esd_dir, pattern = "\\.esd$")

combos <- data.frame(
  celltype = sapply(strsplit(esd_files, "_"), `[`, 1),
  gene = sub("\\.esd$", "", sapply(strsplit(esd_files, "_"), `[`, 2))
)

# ---- 4. Run coloc across every combination ----
all_results <- list()

for (i in seq_len(nrow(combos))) {
  g <- combos$gene[i]
  ct <- combos$celltype[i]
  cat("Running:", g, ct, "\n")
  
  res <- tryCatch(
    run_coloc(g, ct, esd_dir = esd_dir, gwas_dt = gwas_dt),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(res)) {
    all_results[[paste(g, ct, sep = "_")]] <- res
  }
}

# ---- 5. Build and save the summary table ----
summary_table <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    gene = x$gene,
    celltype = x$celltype,
    nsnps = x$summary["nsnps"],
    PP.H4 = x$summary["PP.H4.abf"],
    PP.H3 = x$summary["PP.H3.abf"]
  )
}))

write.csv(summary_table, "results/coloc_results.csv", row.names = FALSE)
cat("\nDone. Results saved to results/coloc_results.csv\n")