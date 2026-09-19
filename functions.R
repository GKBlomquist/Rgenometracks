# ---- imports ----
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(patchwork)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(ggseqlogo)

# ---- set up ----
set_theme(
  theme_classic()+
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12),
      legend.text = element_text(size = 12, color = "black"), legend.title = element_text(size = 12),
      axis.line = element_line(linewidth = 0.4), panel.border = element_rect(linewidth = 0.2, color = "black")
    )
)
# ---- Utility ----
split_region_string <- function(region_string) {
  parts <- str_split(region_string, pattern = ":|-")[[1]]
  chr <- parts[1]
  start <- as.numeric(parts[2])
  end <- as.numeric(parts[3])
  return (list(chr = chr, start = start, end = end))
}


# ---- BigWig -----
fetch_bw <- function(bw_path, chr, start, end) {
  interval <- GRanges(chr, IRanges(start, end))
  bw <- as.data.frame(import.bw(bw_path, which = interval))
  
  return(bw)
}

get_bw_plot <- function(bw_df, start, end, x_axis = FALSE, color = "black", ylabel = "Read Count", xlabel = NULL) {
  plot <- ggplot(bw_df, aes(xmin = start, xmax = end + 1, ymin = 0, ymax = score))+
    geom_rect(fill = color, color = color)+
    scale_y_continuous(expand = c(0,0), name = ylabel)+
    scale_x_continuous(expand = c(0,0), limits = c(start, end))
  
  if (!x_axis) plot = plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  
  return (plot)
}

bigwig_track <- function(file, region, color = "black", ylabel = "Read Count", xlabel = NULL, x_axis = TRUE) {
  region_split <- split_region_string(region)
  bw_df <- fetch_bw(file, region_split[["chr"]], region_split[["start"]], region_split[["end"]])
  plot <- get_bw_plot(bw_df, region_split[["start"]], region_split[["end"]], x_axis = x_axis, color = color, ylabel = ylabel, xlabel = xlabel)
  
  return (plot)
}

# ---- Genes ----
fetch_genes <- function(chr, start, end) {
  interval <- GRanges(chr, IRanges(start, end))
  genedb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  seqlevels(genedb) <- chr
  
  # get exons
  exons <- unlist(exonsBy(genedb, by = "tx"))
  exons$tx_id <- names(exons)
  exons_df <- exons %>% 
    subsetByOverlaps(interval) %>% 
    as.data.frame(row.names = NULL)
  
  # get introns
  introns <- unlist(intronsByTranscript(genedb))
  introns$tx_id <- names(introns)
  introns_df <- introns %>%
    subsetByOverlaps(interval) %>%
    as.data.frame(row.names = NULL)
  
  # group by transcript
  tx_bounds <- exons_df %>%
    group_by(tx_id) %>% 
    summarise(tx_start = min(start), tx_end = max(end)) %>%
    arrange(tx_start)
  
  tx_iranges <- IRanges(start = tx_bounds$tx_start, end = tx_bounds$tx_end)
  tx_bounds$y_level <- disjointBins(tx_iranges)
  
  exons_df <- left_join(exons_df, tx_bounds[, c("tx_id", "y_level")], by = "tx_id")
  introns_df <- left_join(introns_df, tx_bounds[, c("tx_id", "y_level")], by = "tx_id")
  
  return (list(exons = exons_df, introns = introns_df, genes = tx_bounds))
}

annotate_genes <- function(genes) {
  genedb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  tx_to_gene <- select(genedb, 
                       keys = as.character(genes$tx_id), 
                       keytype = "TXID", 
                       columns = "GENEID")
  
  gene_names <- select(org.Hs.eg.db, 
                       keys = tx_to_gene$GENEID, 
                       keytype = "ENTREZID", 
                       columns = "SYMBOL")
  
  genes$symbol <- gene_names$SYMBOL
  
  return(genes)
}

get_genes_plot <- function(introns, exons, start, end, genes = NULL, x_axis = FALSE, transcript = NULL, color = color) {
  if (!is.null(transcript)) {
    introns <- filter(introns, y_level == transcript)
    exons <- filter(exons, y_level == transcript)
    
    if (!is.null(genes)) genes <- filter(genes,  y_level == transcript)
  }
  
  plot <- ggplot() +
    geom_segment(data = introns, 
                 aes(x = start, xend = end, y = y_level, yend = y_level), 
                 color = color, linewidth = 0.5) +
    
    geom_rect(data = exons, 
              aes(xmin = start, xmax = end, 
                  ymin = y_level - 0.2, ymax = y_level + 0.2), 
              fill = color, color = color)+
    
    labs(x = NULL, y = NULL)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
    scale_x_continuous(limits = c(start, end), expand = c(0,0))+
    scale_y_continuous(expand = c(0.05, 0.25))
    
  if (!is.null(genes)) {
    genes <- left_join(genes, exons, by = "y_level") %>% 
      mutate(labelx = ifelse(strand == "+", tx_start - 50, tx_end + 50))
    plot <- plot + geom_text(data = genes, 
                             aes(x = labelx, y = y_level, label = symbol), 
                             hjust = 0, vjust = 0.5, size = 3.5, fontface = "italic")
  
    if (!x_axis) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  }
  return (plot)
}

genes_track <- function(region, x_axis = TRUE, transcript_number = NULL, color = "black") {
  region_split <- split_region_string(region)
  gene_data <- fetch_genes(region_split[["chr"]], region_split[["start"]], region_split[["end"]])
  annotations <- annotate_genes(gene_data[["genes"]])
  plot <- get_genes_plot(gene_data[["introns"]], gene_data[["exons"]], region_split[["start"]], region_split[["end"]], genes = annotations, x_axis = x_axis, transcript = transcript_number, color = color)
  return (plot)
}

# ---- Arcs ----
fetch_arcs <- function(arcs_path, chr, start, end) {
  if (!is.data.frame(arcs_path)) {
    arcs <- read.delim(arcs_path, header = FALSE)
  } else arcs <- arcs_path

  colnames(arcs) <- c("chrom", "start_pos", "end_pos", "score")
  
  arcs_filtered <- arcs %>% 
    filter(chrom == chr & start_pos <= end & end_pos >= start) %>% 
    rowwise() %>% 
    mutate(midpoint = mean(c(start_pos, end_pos))) %>% 
    ungroup() %>%
    mutate(group = row_number()) %>% 
    pivot_longer(cols = c(start_pos, midpoint, end_pos), values_to = "position", names_to = "point")
}

get_arcs_plot <- function(arcs, start, end, inverted = TRUE, x_axis = FALSE, arc_color = "score", color_function = colorRampPalette(c("black", "red")), line_width = 0.5) {
  height = ifelse(inverted, -1, 1)
  
  arcs <- arcs %>% 
    mutate(y = case_when(
      point %in% c("start_pos", "end_pos") ~ 0,
      point == "midpoint" ~ height
    ))
  
  if (arc_color == "score") {
  plot <- ggplot(arcs, aes(x = position, y = y, group = group, color = score))+
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), linewidth = line_width)+
    scale_x_continuous(limits = c(start, end), expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    scale_color_continuous(palette = color_function)+
    labs(x = NULL, y = NULL, color = NULL)+
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      legend.key.height = unit(1, "null"), legend.key.width = unit(0.1, "cm"), legend.frame = element_rect(color = "black")
    )
  } else {
    plot <- ggplot(arcs, aes(x = position, y = y, group = group))+
      geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = arc_color, linewidth = line_width)+
      scale_x_continuous(limits = c(start, end), expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      labs(x = NULL, y = NULL, color = NULL)+
      theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.key.height = unit(1, "null"), legend.key.width = unit(0.1, "cm"), legend.frame = element_rect(color = "black")
      )
  }
  
  
  if (!x_axis) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  
  return (plot)
}

arcs_track <- function(arcs, region, inverted = TRUE, x_axis = TRUE, arc_color = "score", color_function = colorRampPalette(c("black", "red")), line_width = 0.5) {
  region_split <- split_region_string(region)
  arcs <- fetch_arcs(arcs, region_split[["chr"]], region_split[["start"]], region_split[["end"]])
  plot <- get_arcs_plot(arcs, region_split[["start"]], region_split[["end"]], inverted = inverted, x_axis = x_axis, arc_color = arc_color, color_function = color_function, line_width = line_width)
  return (plot)
}

# ---- Logo ----
fetch_contributions <- function(contribution_bw_path, chr, start, end) {
  contributions <- fetch_bw(contribution_bw_path, chr, start, end)
  
  contributions$sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(chr, IRanges(contributions$start, contributions$end)))
  contributions <- dplyr::select(contributions, start, score, sequence) %>% 
    mutate(sequence = as.character(sequence)) %>% 
    pivot_wider(names_from = sequence, values_from = score, values_fill = 0) %>% 
    column_to_rownames("start") %>% 
    t()
  
  return(contributions)
}

fetch_annotations_finemo <- function(finemo_hits_path, contributions, region_chr, region_start , region_end, min_importance = 0, y_padding = 0.1) {
  annotations <- read.delim(finemo_hits_path)

  annotations <- annotations %>% 
    filter(chr == region_chr & start >= region_start & end <= region_end & hit_coefficient_global >= min_importance)
  
  padding = y_padding * max(contributions)
  
  segment_y <- c()
  starts <- c()
  ends <- c()
  for (i in 1:nrow(annotations)) {
    start <- annotations[i, "start"] - (region_start - 1)
    end <- annotations[i, "end"] - (region_start - 1)

    segment_y <- append(segment_y, min(contributions[, start:end]) - padding)
    starts <- append(starts, start)
    ends <- append(ends, end)
  }
  annotations$start <- starts
  annotations$end <- ends
  annotations <- annotations %>% 
    add_column(segment_y = segment_y) %>% 
    rowwise() %>% 
    mutate(label_x = mean(start, end))
  
  return (annotations)
}

get_contribution_plot <- function(contributions, start, end, ylabel = "Contribution", xlabel = NULL, annotations = NULL, x_axis = TRUE) {
  # get breaks
  plot_breaks <- round(seq(0, ncol(contributions), length.out = 4)) + c(1, 0 ,0 ,0)
  
  # plot
  plot <- ggseqlogo(contributions, method = "custom", seq_type ="dna", font="helvetica_bold")+
    labs(y = ylabel, x = xlabel)+
    theme_bw()+
    scale_x_continuous(breaks = plot_breaks, labels = start + plot_breaks, expand = c(0,0))+
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12),
      axis.line = element_line(linewidth = 0.4), panel.border = element_rect(linewidth = 0.2, color = "black")
    )
  
  if (!is.null(annotations)) {
    plot <- plot + 
      geom_segment(data = annotations, aes(x = start, xend = end, y = segment_y, yend = segment_y))+
      geom_text(data = annotations, aes(x = label_x, y = segment_y, label = motif_name), vjust = 1.5, hjust = 0)
  }
  
  if (!x_axis) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  return (plot)
}

bigwig_logo_track <- function(bigwig_path, region, annotations_path = NULL, min_importance = 0, y_padding = 0.1, xlabel = NULL, ylabel = "Contribution", x_axis = TRUE) {
  region_split <- split_region_string(region)
  
  contributions <- fetch_contributions(bigwig_path, region_split[["chr"]], region_split[["start"]], region_split[["end"]])
  
  if (!is.null(annotations_path)) {
    annotations <- fetch_annotations_finemo(annotations_path, contributions, region_split[["chr"]], region_split[["start"]], region_split[["end"]], min_importance = min_importance, y_padding = y_padding)
    }
  else annotations <- NULL
  
  plot <- get_contribution_plot(contributions, region_split[["start"]], region_split[["end"]], ylabel = ylabel, xlabel = xlabel, annotations = annotations, x_axis = x_axis)
}

# ---- test ----
region <- "chr4:82350000-82380000"
zoomed_region <- "chr4:82373900-82374400"
arcs <- data.frame(chr = c("chr4", "chr8"), start = c(82373890, 108443272), end = c(82373990, 108443602), score = c(1, 1))

bigwig_track("~/Documents/ISA/pygenometracks/data/obese_cancer_ATAC.bw", region = region)+
  arcs_track(arcs, region, arc_color = "#E64B35FF")+
  genes_track(region, transcript_number = 4)+
  
  bigwig_track("~/Documents/ISA/chromBPNet/outs/obese/tracks/regulon_sites_prediction_chrombpnet_nobias.bw", zoomed_region, ylabel = "Predicted Profile", x_axis = FALSE)+
  bigwig_logo_track("~/Documents/ISA/chromBPNet/outs/obese/tracks/regulon_sites.profile_scores.bw", zoomed_region, 
                    annotations_path = "Documents/ISA/chromBPNet/outs/obese/finemo/regulon_sites/hits_merged/hits_annotated.tsv", min_importance = 0.000005)+
  
  plot_layout(ncol = 1, axes = "collect", heights = c(1, 1, 0.3, 1, 1)) + plot_annotation(title = "ETS-controlled", theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))) & theme(plot.margin = margin(t = 2, r = 20, b = 2, l = 5, unit = "pt"))


region <- "chr8:108440000-108500000"
zoomed_region <- "chr8:108443400-108443650"

bigwig_track("Documents/ISA/pygenometracks/data/obese_cancer_ATAC.bw", region = region)+
arcs_track(arcs, region, arc_color = "#E64B35FF")+
genes_track(region, transcript_number = 2)+
  
bigwig_track("Documents/ISA/chromBPNet/outs/obese/tracks/regulon_sites_prediction_chrombpnet_nobias.bw", zoomed_region, ylabel = "Predicted Profile", x_axis = FALSE)+
bigwig_logo_track("Documents/ISA/chromBPNet/outs/obese/tracks/regulon_sites.profile_scores.bw", zoomed_region, 
                    annotations_path = "Documents/ISA/chromBPNet/outs/obese/finemo/regulon_sites/hits_merged/hits_annotated.tsv", min_importance = 0.000005)+
  
plot_layout(ncol = 1, axes = "collect", heights = c(1, 1, 0.3, 1, 1)) + plot_annotation(title = "NFY-controlled", theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))) & theme(plot.margin = margin(t = 2, r = 20, b = 2, l = 5, unit = "pt"))

