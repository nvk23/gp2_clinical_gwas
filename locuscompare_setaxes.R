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
# The gene region track under the LocusZoom panels is fetched from the Ensembl
# REST API (GRCh38), so it needs network access. If the lookup fails the script
# warns and draws the figure without the track.
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

gene_track    <- TRUE              # draw the gene region track under the LocusZoom panels
gene_biotypes <- "protein_coding"  # NULL for every biotype, or e.g. c("protein_coding", "lncRNA")

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

# ---- gene annotations ----------------------------------------------------------

#' Genes overlapping the plotted window, from the Ensembl REST API (GRCh38), so
#' no local GTF is needed. Returns a zero-row frame (with a warning) if the
#' lookup fails, which drops the gene track rather than killing the figure.
fetch_genes <- function(chr, start, end, biotypes = "protein_coding") {
  empty <- data.frame(name = character(), start = numeric(),
                      end = numeric(), strand = numeric(),
                      stringsAsFactors = FALSE)

  url <- sprintf(
    "https://rest.ensembl.org/overlap/region/human/%s:%.0f-%.0f?feature=gene;content-type=application/json",
    sub("^chr", "", chr), start, end
  )

  res <- tryCatch(httr::GET(url, httr::timeout(60)), error = function(e) e)
  if (inherits(res, "error") || httr::http_error(res)) {
    warning("Could not fetch gene annotations from Ensembl; skipping the gene track.")
    return(empty)
  }

  g <- tryCatch(
    jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"),
                       simplifyDataFrame = TRUE),
    error = function(e) NULL
  )
  if (!is.data.frame(g) || nrow(g) == 0) return(empty)

  if (!is.null(biotypes) && length(biotypes) && "biotype" %in% names(g)) {
    g <- g[g$biotype %in% biotypes, , drop = FALSE]
  }
  if (nrow(g) == 0) return(empty)

  # prevents unnamed genes from cluttering the track
  nm   <- if ("external_name" %in% names(g)) as.character(g$external_name) else rep(NA_character_, nrow(g))
  keep <- !is.na(nm) & nzchar(nm)
  g    <- g[keep, , drop = FALSE]
  nm   <- nm[keep]
  if (nrow(g) == 0) return(empty)

  out <- data.frame(name   = nm,
                    start  = as.numeric(g$start),
                    end    = as.numeric(g$end),
                    strand = as.numeric(g$strand),
                    stringsAsFactors = FALSE)
  out <- out[order(out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Greedy packing: each gene goes on the topmost row where neither its body nor
#' its (approximate) label width collides with what is already there.
layout_genes <- function(genes, xlim, char_frac = 0.012) {
  if (nrow(genes) == 0) return(cbind(genes, row = integer()))

  span <- diff(xlim)
  genes$start <- pmax(genes$start, xlim[1])
  genes$end   <- pmin(genes$end,   xlim[2])

  mid         <- (genes$start + genes$end) / 2
  half_label  <- nchar(genes$name) * char_frac * span / 2
  genes$left  <- pmin(genes$start, mid - half_label)
  genes$right <- pmax(genes$end,   mid + half_label)

  pad       <- span * 0.01
  row_end   <- numeric(0)
  genes$row <- NA_integer_
  for (i in seq_len(nrow(genes))) {
    placed <- FALSE
    for (r in seq_along(row_end)) {
      if (genes$left[i] > row_end[r] + pad) {
        genes$row[i] <- r
        row_end[r]   <- genes$right[i]
        placed       <- TRUE
        break
      }
    }
    if (!placed) {
      row_end      <- c(row_end, genes$right[i])
      genes$row[i] <- length(row_end)
    }
  }
  genes
}

#' Gene bodies as arrows pointing in the direction of transcription, labelled
#' with the gene symbol. Carries the x axis for the stack above it.
make_gene_track <- function(genes, chr, xlim) {
  genes$x    <- ifelse(genes$strand < 0, genes$end,   genes$start)
  genes$xend <- ifelse(genes$strand < 0, genes$start, genes$end)
  n_rows     <- max(genes$row)

  ggplot(genes) +
    geom_segment(
      aes(x = x, xend = xend, y = row, yend = row),
      linewidth = 0.9, colour = "grey25",
      arrow = arrow(length = unit(0.055, "inches"), type = "closed", ends = "last")
    ) +
    geom_text(
      aes(x = (start + end) / 2, y = row - 0.38, label = name),
      size = 2.8, fontface = "italic", colour = "grey15"
    ) +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x / 1e6),
                       limits = xlim, expand = expansion(mult = 0.02)) +
    scale_y_reverse(limits = c(n_rows + 0.5, 0.35), expand = expansion(mult = 0.05)) +
    xlab(paste0("chr", chr, " (Mb)")) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          plot.margin  = unit(c(0, 1, 0.5, 0.5), "lines"))
}

# ---- panels ------------------------------------------------------------------

#' Lead-SNP label. The rsID is drawn last (so no point covers it), carries a
#' white halo so it stays readable over a dense cloud, and is nudged into the
#' headroom above the point with a leader line back to it.
snp_label <- function(nudge_x = 0, nudge_y = 0) {
  ggrepel::geom_text_repel(
    aes(label = label),
    size               = 3.3,
    fontface           = "bold",
    colour             = "black",
    bg.color           = "white",
    bg.r               = 0.18,
    box.padding        = 0.45,
    point.padding      = 0.45,
    nudge_x            = nudge_x,
    nudge_y            = nudge_y,
    min.segment.length = 0,
    segment.colour     = "grey25",
    segment.size       = 0.4,
    force              = 8,
    max.iter           = 20000,
    max.overlaps       = Inf,
    seed               = 42,
    na.rm              = TRUE
  )
}

make_scatterplot <- function(merged, title1, title2, color, shape, size,
                             nudge_x = 0, nudge_y = 0) {
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
    snp_label(nudge_x, nudge_y) +
    theme_classic()
}

make_locuszoom <- function(stats, title, chr, color, shape, size,
                           xlim = NULL, nudge_y = 0) {
  ggplot(stats, aes(x = pos, logp)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(
      data = stats[stats$label != "", ],
      aes(x = pos, logp, fill = rsid, size = rsid, shape = rsid)
    ) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x / 1e6),
                       limits = xlim, expand = expansion(mult = 0.02)) +
    snp_label(nudge_y = nudge_y) +
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
  # The r^2 title sits above the colour stack, centred on it, so it clears the
  # top box instead of straddling its corner.
  p + cowplot::draw_label(parse(text = "r^2"),
                          x     = legend_box$x[1] + 0.025,
                          y     = max(legend_box$y) + 0.05 + 0.015,
                          hjust = 0.5, vjust = 0, size = 10)
}

#' Assemble the three panels on one shared -log10(P) scale.
make_combined_plot <- function(merged, title1, title2, ld, chr, snp = NULL,
                               r2_legend = TRUE, legend_position = "bottomright",
                               gene_track = TRUE, gene_biotypes = "protein_coding") {
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

  # Headroom at the top of every panel for the nudged rsID label. Continuous
  # limits censor out-of-range values, so y_top has to clear the nudged label
  # position - otherwise ggrepel silently drops the label.
  nudge    <- 0.06 * max_logp
  y_top    <- max_logp + 3 * nudge
  y_expand <- expansion(mult = c(0.02, 0.05))

  p1 <- make_scatterplot(merged, title1, title2, color, shape, size,
                         nudge_y = nudge) +
    scale_x_continuous(limits = c(0, y_top), expand = y_expand) +
    scale_y_continuous(limits = c(0, y_top), expand = y_expand)

  # Must come after the scales: ggdraw() output can no longer take scale_*().
  if (isTRUE(r2_legend)) p1 <- add_r2_legend(p1, legend_position)

  # One x range for every stacked panel, so the gene track lines up with the SNPs.
  xlim <- range(merged$pos, na.rm = TRUE)

  panel <- function(logp_col, title) {
    d <- merged[, c("rsid", logp_col, "chr", "pos", "label")]
    colnames(d)[colnames(d) == logp_col] <- "logp"
    make_locuszoom(d, title, chr, color, shape, size, xlim, nudge_y = nudge) +
      scale_y_continuous(limits = c(0, y_top), expand = y_expand)
  }

  no_x <- theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  p2   <- panel("logp1", title1) + no_x
  p3   <- panel("logp2", title2)

  panels      <- list(p2, p3)
  rel_heights <- c(0.8, 1)

  if (isTRUE(gene_track)) {
    genes <- fetch_genes(chr, xlim[1], xlim[2], gene_biotypes)
    if (nrow(genes) > 0) {
      genes <- layout_genes(genes, xlim)
      # The gene track now carries the x axis, so drop it from the panel above.
      panels      <- list(p2, p3 + no_x, make_gene_track(genes, chr, xlim))
      rel_heights <- c(0.8, 1, 0.22 + 0.1 * max(genes$row))
    }
  }

  right <- cowplot::plot_grid(plotlist = panels, align = "v", axis = "lr",
                              ncol = 1, rel_heights = rel_heights)
  cowplot::plot_grid(p1, right, align = "h", rel_widths = c(1.6, 2))
}

#' Full figure for two already-loaded summary-statistics tables.
locuscompare_shared_axes <- function(stats1, stats2, pos_lookup,
                                     title1, title2, snp = NULL,
                                     population = "EUR", window = 500000,
                                     token = ldlink_token,
                                     r2_legend = TRUE,
                                     legend_position = "bottomright",
                                     gene_track = TRUE,
                                     gene_biotypes = "protein_coding") {
  d1 <- read_sumstats(stats1)
  d2 <- read_sumstats(stats2)

  merged <- merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  merged <- merge(merged, pos_lookup, by = "rsid")

  chr <- unique(merged$chr)
  if (length(chr) != 1) stop("There must be one and only one chromosome.")

  snp <- get_lead_snp(merged, snp)
  ld  <- retrieve_LD(snp, population, window, token)
  make_combined_plot(merged, title1, title2, ld, chr, snp,
                     r2_legend = r2_legend, legend_position = legend_position,
                     gene_track = gene_track, gene_biotypes = gene_biotypes)
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
  legend_position = legend_position,
  gene_track      = gene_track,
  gene_biotypes   = gene_biotypes
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir,
                      sprintf("locuscompare_chr%s_%s_%s_%s.png", chr, bp, ref, alt))
ggsave(filename = out_file, plot = p, width = 12, height = 6, dpi = 300)
cat("wrote", out_file, "\n")
