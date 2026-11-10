#!/usr/bin/env Rscript
#
# LocusCompare plot with shared axis limits.
#
# Draws the three-panel locuscompare figure for one locus: risk GWAS against
# progression GWAS as a scatter, plus a LocusZoom panel for each. Unlike the
# stock locuscomparer, all panels share one -log10(P) scale, so the two studies
# can be read against each other directly.
#
# Inputs are the LocusZoom-format summary statistics written by
# extra/locuszoom_prep.ipynb: bgzipped, tab separated, with
# columns chrom, pos, rsid, neg_log_pvalue.
#
# Usage:
#   export LDLINK_TOKEN=...          # free: https://ldlink.nih.gov/?tab=apiaccess
#   export LC_RISK_FILE=...          # risk (incidence) GWAS
#   export LC_PROGRESSION_FILE=...   # progression GWAS
#   export LC_OUT_DIR=...            # where the PNG is written
#   Rscript locuscompare_setaxes.R
#
# The locus and plotting options are set in the block below.

library(locuscomparer)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(LDlinkR)

# ---- locus to plot (hg38) ----
chr      <- "13"
bp       <- 47362415
ref      <- "A"
alt      <- "G"
lead_snp <- "rs55833024"
window   <- 500000            # +/- bp around `bp`
population <- "EUR"           # 1000 Genomes population used for LD

title_risk        <- "RISK"
title_progression <- "SURVIVAL"

r2_legend       <- TRUE            # draw the r^2 colour key on the scatter panel
legend_position <- "bottomright"   # "bottomright", "topright" or "topleft"

# ---- inputs from the environment ----
ldlink_token     <- Sys.getenv("LDLINK_TOKEN")
risk_file        <- Sys.getenv("LC_RISK_FILE")
progression_file <- Sys.getenv("LC_PROGRESSION_FILE")
out_dir          <- Sys.getenv("LC_OUT_DIR", unset = "locuscompare_output")

for (v in c("risk_file", "progression_file")) {
  if (!nzchar(get(v))) {
    stop(sprintf("set LC_%s in your environment (see the header of this script)",
                 toupper(sub("_file$", "_FILE", v))))
  }
}

# ---- data preparation --------------------------------------------------------

#' Read one summary-statistics table down to rsid / pval / logp.
read_sumstats <- function(x, marker_col = "rsid", pval_col = "pval") {
  if (is.character(x)) {
    d <- read.table(x, header = TRUE, stringsAsFactors = FALSE)
    colnames(d)[colnames(d) == marker_col] <- "rsid"
    colnames(d)[colnames(d) == pval_col]   <- "pval"
  } else if (is.data.frame(x)) {
    d <- x
  } else {
    stop("x must be a file path or a data.frame")
  }
  d$logp <- -log10(d$pval)
  d[, c("rsid", "pval", "logp")]
}

#' awk one locus out of a bgzipped LocusZoom file, so only the region is read.
region_cmd <- function(fn, chr, start, end) {
  sprintf("gzip -dc %s | awk -F'\\t' 'NR==1 || ($1==\"%s\" && $2>=%.0f && $2<=%.0f)'",
          shQuote(fn), chr, start, end)
}

#' Read one locus and split it into the p-value table and a position lookup.
read_locus <- function(fn, chr, start, end) {
  d <- fread(cmd = region_cmd(fn, chr, start, end))
  setnames(d, 1, "chrom")
  d <- d[rsid != "." & !is.na(neg_log_pvalue)]
  d[, pval := 10^(-neg_log_pvalue)]
  d <- unique(d, by = "rsid")
  list(
    stats = as.data.frame(d[, .(rsid, pval)]),
    pos   = d[, .(rsid, chr = as.character(chrom), pos = as.numeric(pos))]
  )
}

#' LD (r^2) between the lead SNP and everything around it, from the LDlink API.
retrieve_LD <- function(snp, population, window, token) {
  if (!nzchar(token)) {
    stop("No LDlink token. Get one at https://ldlink.nih.gov/?tab=apiaccess ",
         "and export LDLINK_TOKEN=...")
  }
  # NIH serves the idp.nih.gov certificate on ldlink.nih.gov, so the hostname
  # check fails. Skip that check only; ssl_verifypeer stays on, so the chain is
  # still validated against a trusted CA.
  res <- httr::with_config(
    httr::config(ssl_verifyhost = 0),
    LDlinkR::LDproxy(
      snp          = snp,
      pop          = population,
      r2d          = "r2",
      token        = token,
      genome_build = "grch38",
      win_size     = window
    )
  )
  if (!"RS_Number" %in% colnames(res)) {
    stop(sprintf("LDproxy returned no LD for %s: %s", snp,
                 paste(utils::capture.output(print(utils::head(res, 2))),
                       collapse = " | ")))
  }
  res <- res[res$RS_Number != "." & !is.na(res$R2), ]
  data.frame(SNP_A = snp,
             SNP_B = as.character(res$RS_Number),
             R2    = as.numeric(res$R2),
             stringsAsFactors = FALSE)
}

#' Lead SNP: the one given, or the strongest joint signal across both studies.
get_lead_snp <- function(merged, snp = NULL) {
  if (is.null(snp)) {
    snp <- merged[which.min(merged$pval1 + merged$pval2), "rsid"]
  } else if (!snp %in% merged$rsid) {
    stop(sprintf("%s is not in the intersection of the two inputs.", snp))
  }
  as.character(snp)
}

#' Standard locuscompare r^2 bins; the lead SNP is purple, unknown LD is blue4.
assign_color <- function(rsid, snp, ld) {
  ld <- ld[ld$SNP_A == snp, ]
  ld$color <- as.character(cut(
    ld$R2,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c("blue4", "skyblue", "darkgreen", "orange", "red"),
    include.lowest = TRUE
  ))

  color <- data.frame(rsid, stringsAsFactors = FALSE)
  color <- merge(color, ld[, c("SNP_B", "color")],
                 by.x = "rsid", by.y = "SNP_B", all.x = TRUE)
  color[is.na(color$color), "color"] <- "blue4"

  if (snp %in% color$rsid) {
    color[color$rsid == snp, "color"] <- "purple"
  } else {
    color <- rbind(color, data.frame(rsid = snp, color = "purple"))
  }

  res <- color$color
  names(res) <- color$rsid
  res
}

add_label <- function(merged, snp) {
  merged$label <- ifelse(merged$rsid %in% snp, merged$rsid, "")
  merged
}

# ---- panels ------------------------------------------------------------------

make_scatterplot <- function(merged, title1, title2, color, shape, size) {
  ggplot(merged, aes(logp1, logp2)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(
      data = merged[merged$label != "", ],
      aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)
    ) +
    xlab(bquote(.(title1) ~ -log[10] * "(P)")) +
    ylab(bquote(.(title2) ~ -log[10] * "(P)")) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf) +
    theme_classic()
}

make_locuszoom <- function(stats, title, chr, color, shape, size) {
  ggplot(stats, aes(x = pos, logp)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(
      data = stats[stats$label != "", ],
      aes(x = pos, logp, fill = rsid, size = rsid, shape = rsid)
    ) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x / 1e6)) +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf) +
    xlab(paste0("chr", chr, " (Mb)")) +
    ylab(bquote(.(title) ~ -log[10] * "(P)")) +
    theme_classic() +
    theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"))
}

add_r2_legend <- function(p, position = c("bottomright", "topright", "topleft")) {
  position <- match.arg(position)
  legend_box <- if (position == "bottomright") {
    data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
  } else if (position == "topright") {
    data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
  } else {
    data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
  }

  p <- cowplot::ggdraw(p) +
    geom_rect(
      data = legend_box,
      aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
      color = "black",
      fill  = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))
    )
  for (i in seq_len(4)) {
    p <- p + cowplot::draw_label(c("0.8", "0.6", "0.4", "0.2")[i],
                                 x = legend_box$x[i] + 0.05,
                                 y = legend_box$y[i], hjust = -0.3, size = 10)
  }
  p + cowplot::draw_label(parse(text = "r^2"),
                          x = legend_box$x[1] + 0.05, y = legend_box$y[1],
                          vjust = -2, size = 10)
}

#' Assemble the three panels on one shared -log10(P) scale.
make_combined_plot <- function(merged, title1, title2, ld, chr, snp = NULL,
                               r2_legend = TRUE, legend_position = "bottomright") {
  snp   <- get_lead_snp(merged, snp)
  color <- assign_color(merged$rsid, snp, ld)

  shape <- ifelse(merged$rsid == snp, 23, 21)
  names(shape) <- merged$rsid

  size <- ifelse(merged$rsid == snp, 3, 2)
  names(size) <- merged$rsid

  merged <- add_label(merged, snp)

  max_logp <- max(c(merged$logp1, merged$logp2), na.rm = TRUE)
  if (!is.finite(max_logp) || max_logp <= 0) max_logp <- 1
  max_logp <- ceiling(max_logp)

  p1 <- make_scatterplot(merged, title1, title2, color, shape, size) +
    scale_x_continuous(limits = c(0, max_logp)) +
    scale_y_continuous(limits = c(0, max_logp))

  # Must come after the scales: ggdraw() output can no longer take scale_*().
  if (isTRUE(r2_legend)) p1 <- add_r2_legend(p1, legend_position)

  panel <- function(logp_col, title) {
    d <- merged[, c("rsid", logp_col, "chr", "pos", "label")]
    colnames(d)[colnames(d) == logp_col] <- "logp"
    make_locuszoom(d, title, chr, color, shape, size) +
      scale_y_continuous(limits = c(0, max_logp))
  }

  p2 <- panel("logp1", title1) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  p3 <- panel("logp2", title2)

  right <- cowplot::plot_grid(p2, p3, align = "v", nrow = 2,
                              rel_heights = c(0.8, 1))
  cowplot::plot_grid(p1, right, align = "h", rel_widths = c(1.6, 2))
}

#' Full figure for two already-loaded summary-statistics tables.
locuscompare_shared_axes <- function(stats1, stats2, pos_lookup,
                                     title1, title2, snp = NULL,
                                     population = "EUR", window = 500000,
                                     token = ldlink_token,
                                     r2_legend = TRUE,
                                     legend_position = "bottomright") {
  d1 <- read_sumstats(stats1)
  d2 <- read_sumstats(stats2)

  merged <- merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  merged <- merge(merged, pos_lookup, by = "rsid")

  chr <- unique(merged$chr)
  if (length(chr) != 1) stop("There must be one and only one chromosome.")

  snp <- get_lead_snp(merged, snp)
  ld  <- retrieve_LD(snp, population, window, token)
  make_combined_plot(merged, title1, title2, ld, chr, snp,
                     r2_legend = r2_legend, legend_position = legend_position)
}

# ---- run ---------------------------------------------------------------------

start <- max(0, bp - window)
end   <- bp + window

risk <- read_locus(risk_file, chr, start, end)
prog <- read_locus(progression_file, chr, start, end)

pos_lookup <- as.data.frame(unique(rbind(risk$pos, prog$pos), by = "rsid"))

p <- locuscompare_shared_axes(
  stats1          = risk$stats,
  stats2          = prog$stats,
  pos_lookup      = pos_lookup,
  title1          = title_risk,
  title2          = title_progression,
  snp             = lead_snp,
  population      = population,
  window          = window,
  r2_legend       = r2_legend,
  legend_position = legend_position
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir,
                      sprintf("locuscompare_chr%s_%s_%s_%s.png", chr, bp, ref, alt))
ggsave(filename = out_file, plot = p, width = 12, height = 6, dpi = 300)
cat("wrote", out_file, "\n")
