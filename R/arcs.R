# ---- Utility ----
split_region_string <- function(region_string) {
  parts <- str_split(region_string, pattern = ":|-")[[1]]
  chr <- parts[1]
  start <- as.numeric(parts[2])
  end <- as.numeric(parts[3])
  return (list(chr = chr, start = start, end = end))
}


# ---- IO ----
fetch_arcs <- function(arcs_path, chr, start, end) {
  if (!is.data.frame(arcs_path)) {
    arcs <- read.delim(arcs_path, header = FALSE)
  } else arcs <- arcs_path

  colnames(arcs) <- c("chrom", "start_pos", "end_pos", "score")

  arcs_filtered <- arcs %>%
    dplyr::filter(chrom == chr & start_pos <= end & end_pos >= start) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(midpoint = mean(c(start_pos, end_pos))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = c(start_pos, midpoint, end_pos), values_to = "position", names_to = "point")

  return(arcs_filtered)
}


# ---- Plotting ----
get_arcs_plot <- function(arcs, start, end, inverted = TRUE, x_axis = FALSE, arc_color = "score", color_function = grDevices::colorRampPalette(c("black", "red")), line_width = 0.5) {
  height = ifelse(inverted, -1, 1)

  arcs <- arcs %>%
    dplyr::mutate(y = dplyr::case_when(
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


# ---- Exposed functions ----
#' @title title
#' @description
#' A short description...
#' @param arcs, character | data.frame: Path to bed file with the arcs. Should have the following format: chr, start, end, score. Alternatively a data.frame can be supplied directly with the appropriate values. Column names are not considered. Thus, the columns should be in the right order.
#' @param region, character: string specifying region to plot.
#' @param inverted, bool: whether to invert the arcs.
#' @param x_axis, bool: whether to draw the x-axis.
#' @param arc_color, character: color of the arcs. If 'score' the arcs will be colored according to the score column.
#' @param color_function, function: color palette function used to map scores to colors.
#' @param line_width, numeric: width of the arc lines.
#' @returns A ggplot2 plot.
#' @import ggplot2
#' @export
arcs_track <- function(arcs, region, inverted = TRUE, x_axis = TRUE, arc_color = "score", color_function = grDevices::colorRampPalette(c("black", "red")), line_width = 0.5) {
  region_split <- split_region_string(region)
  arcs_data <- fetch_arcs(arcs, region_split[["chr"]], region_split[["start"]], region_split[["end"]])
  plot <- get_arcs_plot(arcs_data, region_split[["start"]], region_split[["end"]], inverted = inverted, x_axis = x_axis, arc_color = arc_color, color_function = color_function, line_width = line_width)
  return (plot)
}
