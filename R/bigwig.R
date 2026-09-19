# ---- Utility ----
split_region_string <- function(region_string) {
  parts <- stringr::str_split(region_string, pattern = ":|-")[[1]]
  chr <- parts[1]
  start <- as.numeric(parts[2])
  end <- as.numeric(parts[3])
  return (list(chr = chr, start = start, end = end))
}


# ---- IO ----
fetch_bw <- function(bw_path, chr, start, end) {
  interval <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  bw <- as.data.frame(rtracklayer::import.bw(bw_path, which = interval))

  return(bw)
}


fetch_contributions <- function(contribution_bw_path, chr, start, end, genome) {
  contributions <- fetch_bw(contribution_bw_path, chr, start, end)

  contributions$sequence <- Biostrings::getSeq(genome, GenomicRanges::GRanges(chr, IRanges::IRanges(contributions$start, contributions$end)))
  contributions <- dplyr::select(contributions, start, score, sequence) %>%
    dplyr::mutate(sequence = as.character(sequence)) %>%
    tidyr::pivot_wider(names_from = sequence, values_from = score, values_fill = 0) %>%
    tibble::column_to_rownames("start") %>%
    t()

  return(contributions)
}


fetch_annotations_finemo <- function(finemo_hits_path, contributions, region_chr, region_start , region_end, min_importance = 0, y_padding = 0.1) {
  annotations <- read.delim(finemo_hits_path)

  annotations <- annotations %>%
    dplyr::filter(chr == region_chr & start >= region_start & end <= region_end & hit_coefficient_global >= min_importance)

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
    tibble::add_column(segment_y = segment_y) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(label_x = mean(c(start, end)))

  return (annotations)
}


# ---- Plotting ----
get_bw_plot <- function(bw_df, start, end, x_axis = FALSE, color = "black", ylabel = "Read Count", xlabel = NULL, y_max = NULL) {
  if (is.null(y_max)) {
    y_max <- max(bw_df$score)
  }

  plot <- ggplot(bw_df, aes(xmin = start, xmax = end + 1, ymin = 0, ymax = score))+
    geom_rect(fill = color, color = color)+
    scale_y_continuous(expand = c(0,0), name = ylabel, limits = c(0, y_max))+
    scale_x_continuous(expand = c(0,0), limits = c(start, end))

  if (!x_axis) plot = plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

  return (plot)
}


get_contribution_plot <- function(contributions, start, end, ylabel = "Contribution", xlabel = NULL, annotations = NULL, x_axis = TRUE, y_max = NULL) {
  # get breaks
  plot_breaks <- round(seq(0, ncol(contributions), length.out = 4)) + c(1, 0 ,0 ,0)

  # get y-limit
  if (is.null(y_max)) {
    y_max <- max(contributions)
  }

  # plot
  plot <- ggseqlogo::ggseqlogo(contributions, method = "custom", seq_type ="dna", font="helvetica_bold")+
    labs(y = ylabel, x = xlabel)+
    theme_bw()+
    scale_x_continuous(breaks = plot_breaks, labels = start + plot_breaks, expand = c(0,0))+
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12),
      axis.line = element_line(linewidth = 0.4), panel.border = element_rect(linewidth = 0.2, color = "black")
    )+
    ylim(c(0, y_max))

  if (!is.null(annotations)) {
    plot <- plot +
      geom_segment(data = annotations, aes(x = start, xend = end, y = segment_y, yend = segment_y))+
      geom_text(data = annotations, aes(x = label_x, y = segment_y, label = motif_name), vjust = 1.5, hjust = 0)
  }

  if (!x_axis) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  return (plot)
}


# ---- exposed functions ----
#' @title Plot bigwig track at specified chromosomal region
#' @description
#' Plot specified region of bigwig file using ggplot2.
#' @param file, character: path to BigWig file.
#' @param region, character: string specifying region to plot (accepts chr:start-end and chr-start-end formats).
#' @param color, character: fill color of plot.
#' @param ylabel, character| NULL: y-axis label.
#' @param xlabel, character | NULL: x-axis label.
#' @param x_axis, bool: Whether to draw the x-axis.
#' @param y_max, numeric: Upper limit of the y-axis.
#' @returns A ggplot 2 plot.
#' @import ggplot2
#' @export
bigwig_track <- function(file, region, color = "black", ylabel = "Read Count", xlabel = NULL, x_axis = TRUE, y_max = NULL) {
  region_split <- split_region_string(region)
  bw_df <- fetch_bw(file, region_split[["chr"]], region_split[["start"]], region_split[["end"]])
  plot <- get_bw_plot(bw_df, region_split[["start"]], region_split[["end"]], x_axis = x_axis, color = color, ylabel = ylabel, xlabel = xlabel, y_max = y_max)

  return (plot)
}


#' @title Plot attribution score bigwig logo track at specified chromosomal region
#' @description
#' Plots bigwig scores at chromosomal region as a logo track, where the height of each letter corresponds to the score at the respective coordinate. Intended for deep learning model attribution scores.
#' @param file, character: path to BigWig file.
#' @param region, character: string specifying region to plot (accepts chr:start-end and chr-start-end formats).
#' @param BSgenome, BSgenome: A BSgenome object for the reference genome used. For instance pass BSgenome.Hsapiens.UCSC.hg38 for hg38 reference genome.
#' @param annotations_path, character | NULL: Path to annotation file. If provided, bars and names will be drawn under respective regions. Intended for finemo motif-search output, but will accept any tsv file with columns region_chr, region_start, region_end, hit_coefficient_global, and motif_name.
#' @param min_importance, numeric: Minimum importance-score cutoff for annotations. Any annotation under this threshold will not be drawn.
#' @param y_padding, numeric: Spacing between logo letters and annotations.
#' @param ylabel, character| NULL: y-axis label.
#' @param xlabel, character | NULL: x-axis label.
#' @param x_axis, bool: Whether to draw the x_axis.
#' @param y_max, numeric: Upper limit of the y-axis.
#' @returns A ggplot 2 plot.
#' @import ggplot2
#' @export
bigwig_logo_track <- function(file, region, BSgenome, annotations_path = NULL, min_importance = 0, y_padding = 0.1, xlabel = NULL, ylabel = "Contribution", x_axis = TRUE, y_max = NULL) {
  region_split <- split_region_string(region)

  contributions <- fetch_contributions(file, region_split[["chr"]], region_split[["start"]], region_split[["end"]], BSgenome)

  if (!is.null(annotations_path)) {
    annotations <- fetch_annotations_finemo(annotations_path, contributions, region_split[["chr"]], region_split[["start"]], region_split[["end"]], min_importance = min_importance, y_padding = y_padding)
  }
  else annotations <- NULL

  plot <- get_contribution_plot(contributions, region_split[["start"]], region_split[["end"]], ylabel = ylabel, xlabel = xlabel, annotations = annotations, x_axis = x_axis, y_max = y_max)
  return(plot)
}
