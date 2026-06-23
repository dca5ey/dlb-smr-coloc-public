library(coloc)
library(ggplot2)

run_coloc <- function(gene, celltype, esd_dir, gwas_dt,
                      eqtl_N = 424, gwas_case_n = 2981, gwas_control_n = 4390) {
  
  eqtl_file <- file.path(esd_dir, paste0(celltype, "_", gene, ".esd"))
  eqtl_data <- read.delim(eqtl_file, header = TRUE, sep = "\t")
  
  gwas_subset <- gwas_dt[SNP %in% eqtl_data$SNP]
  merged <- merge(eqtl_data, gwas_subset, by = "SNP", suffixes = c(".eqtl", ".gwas"))
  
  eqtl_input <- list(
    beta = merged$Beta, varbeta = merged$se.eqtl^2,
    MAF = pmin(merged$Freq, 1 - merged$Freq),
    N = eqtl_N, type = "quant",
    snp = merged$SNP, position = merged$Bp
  )
  
  gwas_input <- list(
    beta = merged$b, varbeta = merged$se.gwas^2,
    MAF = pmin(merged$freq, 1 - merged$freq),
    N = merged$n, s = gwas_case_n / (gwas_case_n + gwas_control_n),
    type = "cc", snp = merged$SNP, position = merged$Bp
  )
  
  result <- coloc.abf(dataset1 = eqtl_input, dataset2 = gwas_input)
  list(gene = gene, celltype = celltype, nsnps = nrow(merged), summary = result$summary)
}

make_locus_plot <- function(gene, celltype) {
  eqtl_file <- paste0("/Volumes/RESEARCH/DLB_eQTL/data/eqtl/besd_input/", celltype, "_", gene, ".esd")
  eqtl_data <- read.delim(eqtl_file, header = TRUE, sep = "\t")
  gwas_subset <- gwas_dt[SNP %in% eqtl_data$SNP]
  merged <- merge(eqtl_data, gwas_subset, by = "SNP", suffixes = c(".eqtl", ".gwas"))
  
  plot_data <- rbind(
    data.frame(Bp = merged$Bp, log10p = -log10(merged$p.eqtl), track = "eQTL"),
    data.frame(Bp = merged$Bp, log10p = -log10(merged$p.gwas), track = "GWAS")
  )
  
  ggplot(plot_data, aes(x = Bp / 1e6, y = log10p)) +
    geom_point(alpha = 0.6, size = 1, color = "steelblue") +
    facet_wrap(~track, ncol = 1, scales = "free_y") +
    labs(title = paste(gene, "in", celltype), x = "Position (Mb)", y = "-log10(p)") +
    theme_minimal()
}