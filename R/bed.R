# ---- Utility ----
split_region_string <- function(region_string) {
  parts <- stringr::str_split(region_string, pattern = ":|-")[[1]]
  chr <- parts[1]
  start <- as.numeric(parts[2])
  end <- as.numeric(parts[3])
  return (list(chr = chr, start = start, end = end))
}


# ---- IO ----
fetch_bed <- function(bed_path, chr, start, end) {
  interval <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  bed <- as.data.frame(rtracklayer::import(bed_path, which = interval))

  return(bed)
}


# ---- Plotting ----
get_bed_plot <- function(beds_df, start, end, x_axis = TRUE, color = NULL, xlabel = NULL, legend = TRUE, linewidth = 2, ylabel = NULL, y_axis = TUE) {
  # base plot
  plot <- ggplot(beds_df, aes(x = start, xend = end, y = sample, yend = sample)) +
    scale_x_continuous(expand = c(0,0), limits = c(start, end), name = xlabel) +
    scale_y_discrete(name = ylabel)

  # color options
    if (is.character(color) | length(unique(beds_df$sample)) == 1) {

      if (is.null(color)) color = "black"
      plot <- plot + geom_segment(color = color, linewidth = linewidth)

    } else {
      plot <- plot + geom_segment(aes(color = sample), show.legend = legend, linewidth = linewidth)

      if (!is.null(color)) plot <- plot + scale_color_discrete(palette = color)
    }

  if (!x_axis) {
    plot <- plot + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank()
    )
  }

  if (!y_axis) {
    plot <- plot + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
  }

  return (plot)
}

# ---- Exported functions ----
#' @title Plot one or multiple bed tracks in same region.
#' @description
#' Plots a region over multiple bed tracks.
#' @param bed_paths, named list | character: Either a path to a single bed file or a named list of bed files, where each entry in the list is the path to the file, and the name being the sample name (e.g., list('sample1' = 'path/to/sample1.bed', 'sample2' = 'path/to/sample2.bed').
#' @param region, character: string specifying region to plot (accepts chr:start-end and chr-start-end formats).
#' @param colors, NULL | vector | function: Either a character vector of colors to use (must be at least same length as number of files) or a color function determining how samples are colored. If NULL, will default to ggplot2 default palette.
#' @param ylabel, character | NULL: y-axis label.
#' @param xlabel, character | NULL: x-axis label.
#' @param x_axis, bool: Whether to draw the x_axis.
#' @param y_axis, bool: Whether to draw the y_axis.
#' @param legend, bool: Whether to draw color legend.
#' @param linewidth, numeric: Width of the lines in the bed plot.
#' @returns A ggplot 2 plot.
#' @import ggplot2
#' @export
bed_track <- function(bed_paths, region, x_axis = TRUE, color = NULL, xlabel = NULL, legend = TRUE, linewidth = 2, ylabel = NULL, y_axis = TRUE) {
  region_split <- split_region_string(region)

  if (is.character(bed_paths)) {
    bed_df <- fetch_bed(bed_paths, region_split[["chr"]], region_split[["start"]], region_split[["end"]])
    bed_df$sample = ""
  } else {

    bed_list <- list()
    for (sample in names(bed_paths)) {
      bed_list[[sample]] <- fetch_bed(bed_paths[[sample]], region_split[["chr"]], region_split[["start"]], region_split[["end"]])
    }

    bed_df <- dplyr::bind_rows(bed_list, .id = "sample")
  }

  plot <- get_bed_plot(bed_df, region_split[["start"]], region_split[["end"]], x_axis = x_axis, color = color, xlabel = xlabel, legend = legend, linewidth = linewidth, ylabel = ylabel, y_axis = y_axis)
  return (plot)
}
