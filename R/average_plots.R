fetch_one_average_profile <- function(bw_path, peaks_path, padding, bin_size) {
  # ---- read peaks ----
  peaks_df <- read.delim(peaks_path, header = FALSE)[, 1:3]
  colnames(peaks_df) <- c("chr", "start", "end")

  peaks <- GenomicRanges::GRanges(
    seqnames = peaks_df$chr,
    ranges = IRanges::IRanges(start = peaks_df$start, end = peaks_df$end)
  )

  # ---- import bigwig and sync seqinfo ----
  bw_file <- rtracklayer::BigWigFile(bw_path)
  bw_seqinfo <- GenomeInfoDb::seqinfo(bw_file)

  GenomeInfoDb::seqlevels(peaks, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(bw_seqinfo)
  GenomeInfoDb::seqinfo(peaks) <- bw_seqinfo

  peak_centers <- GenomicRanges::start(peaks) + floor(GenomicRanges::width(peaks) / 2)

  # ---- extend ranges using padding ----
  window_starts <- peak_centers - padding
  window_ends <- peak_centers + padding

  windows <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(peaks),
    ranges = IRanges::IRanges(start = window_starts, end = window_ends),
    seqinfo = bw_seqinfo
  )

  windows <- GenomicRanges::trim(windows)

  # ---- import bigwig for reduced windows ----
  bw <- rtracklayer::import.bw(bw_file, which = GenomicRanges::reduce(windows), as = "RleList")

  # ---- tile and get average ----
  bw_binned <- unlist(GenomicRanges::tile(windows, width = bin_size))
  GenomeInfoDb::seqinfo(bw_binned) <- bw_seqinfo
  average_bins <- GenomicRanges::binnedAverage(bw_binned, bw, "score", na.rm = TRUE)

  n_peaks <- length(peaks)
  score_matrix <- matrix(average_bins$score, nrow = n_peaks, byrow = TRUE)
  rownames(score_matrix) <- 1:nrow(score_matrix)
  ranked_order <- names(rev(sort(rowMeans(score_matrix))))
  return(score_matrix[ranked_order, ])
}


#' @title Fetch average bigwig profiles over peak regions
#' @description
#' Returns bigwig data centered around the input peaks at a specified window and bin size. Output data depends on the parameter plot_type. If plot_typew is 'line' the function will return a data frame with 'position', 'score' and 'sample_name' if multiple samples are inputted.
#' @param bw_path, named list | character: A named list of bigwig files, where each entry in the list is the path to the file, and the name being the sample name (e.g., list('sample1' = 'path/to/sample1.bw', 'sample2' = 'path/to/sample2.bw'). or alternatively path to a single bigwig fdile
#' @param peaks_path, character: path to bed file with regions).
#' @param padding, integer: number of base pairs to extend from the center of the provided peaks.
#' @param bin_size, integer: number of bins to use.
#' @param plot_type, character: desired plot type. Must be either 'line' or 'heatmap'. Will influence the structure of the output data. If 'line', bigwig scores will be averaged across the regions for each bin.
#' @returns Either a data frame of the average bigwig score per sample per position, or a named list of matrices for each sample (or single matrix if only one sample) of bin-wise scores (columns) for each region (rows).
#' @export
fetch_average_profile <- function(bw_path, peaks_path, padding = 1000, bin_size = 10, plot_type = "line") {
  if (is.character(bw_path)) {
    score_matrix <- fetch_one_average_profile(bw_path, peaks_path, padding, bin_size)

    if (plot_type == "line") {
      out_data <- data.frame(score = colMeans(score_matrix),
                             position = seq(-padding, padding, by = bin_size))
    } else if (plot_type == "heatmap") {
      out_data <- score_matrix
    } else {
      stop("plot_type must be either 'line' or 'heatmap'")
    }

  } else if (is.list(bw_path)) {

    if (plot_type == "line") {
      out_data <- dplyr::bind_rows(purrr::imap(bw_path, function(path, sample_name) {
        score_matrix <- fetch_one_average_profile(path, peaks_path, padding, bin_size)

        data_frame <- data.frame(
          score = colMeans(score_matrix),
          position = seq(-padding, padding, by = bin_size)
        )
      }, .progress = TRUE), .id = "sample_name")

    } else if (plot_type == "heatmap") {
      out_data <- purrr::imap(bw_path, function(path, sample_name) {
        fetch_one_average_profile(path, peaks_path, padding, bin_size)
      }, .progress = TRUE)
    }

  } else {
    stop("bw_path must be either character or named list of character")
  }

  return (out_data)
}
