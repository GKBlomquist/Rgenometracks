# ---- IO ----
# this function is currently a disgusting mess
fetch_genes <- function(chr, chrom_start, chrom_end, txdb, orgdb, collapse = FALSE, fully_in_view = FALSE) {
  interval <- GenomicRanges::GRanges(chr, IRanges::IRanges(chrom_start, chrom_end))
  GenomeInfoDb::seqlevels(txdb) <- chr

  visible_tx <- GenomicFeatures::transcriptsByOverlaps(txdb, interval)
  tx_df <- as.data.frame(visible_tx)

  if (nrow(tx_df) == 0) {
    return(list(exons = data.frame(), introns = data.frame(), labels = data.frame()))
  }

  # Filter for fully in view genes
  if (fully_in_view) {
    tx_df <- tx_df %>% dplyr::filter(start >= chrom_start & end <= chrom_end)
  }

  if (nrow(tx_df) == 0) {
    return(list(exons = data.frame(), introns = data.frame(), labels = data.frame()))
  }

  # Annotate transcripts with Gene Symbols
  tx_ids <- as.character(tx_df$tx_id)

  suppressMessages({
    tx_to_gene <- AnnotationDbi::select(txdb, keys = tx_ids, keytype = "TXID", columns = "GENEID")
    tx_to_gene <- tx_to_gene[!duplicated(tx_to_gene$TXID) & !is.na(tx_to_gene$GENEID), ]

    if (nrow(tx_to_gene) > 0) {
      gene_names <- AnnotationDbi::select(orgdb, keys = unique(tx_to_gene$GENEID), keytype = "ENTREZID", columns = "SYMBOL")
      gene_names <- gene_names[!duplicated(gene_names$ENTREZID), ]
      tx_map <- dplyr::left_join(tx_to_gene, gene_names, by = c("GENEID" = "ENTREZID"))
    } else {
      tx_map <- data.frame(TXID = character(), SYMBOL = character())
    }
  })

  tx_df$TXID <- as.character(tx_df$tx_id)
  tx_df <- dplyr::left_join(tx_df, tx_map, by = "TXID")

  # Remove any transcripts that dont have a known gene symbol
  tx_df <- tx_df[!is.na(tx_df$SYMBOL), ]

  # return empty data if no annotated transcripts were found
  if (nrow(tx_df) == 0) {
    return(list(exons = data.frame(), introns = data.frame(), labels = data.frame()))
  }

  # Fetch exons and introns
  exons_list <- GenomicFeatures::exonsBy(txdb, by = "tx")[tx_ids]
  introns_list <- GenomicFeatures::intronsByTranscript(txdb)[tx_ids]

  # Handle grouping based on collapse argument
  if (collapse) {
    # Attach symbols to the GRanges object before reducing
    exons_gr <- unlist(exons_list)
    tx_to_sym <- stats::setNames(tx_df$SYMBOL, tx_df$TXID)
    exons_gr$SYMBOL <- tx_to_sym[names(exons_gr)]

    # Split by symbol and merge overlapping exons
    exons_split <- S4Vectors::split(exons_gr, exons_gr$SYMBOL)
    exons_reduced <- GenomicRanges::reduce(exons_split)

    exons_df <- as.data.frame(exons_reduced) %>% dplyr::rename(SYMBOL = group_name)

    # calculate introns as the gaps between merged exons using dplyr
    introns_df <- exons_df %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::arrange(start, .by_group = TRUE) %>%
      dplyr::mutate(
        intron_start = end,
        intron_end = dplyr::lead(start)
      ) %>%
      dplyr::filter(!is.na(intron_end)) %>%
      dplyr::mutate(
        start = intron_start,
        end = intron_end,
        width = end - start
      ) %>%
      dplyr::select(-intron_start, -intron_end) %>%
      dplyr::ungroup()

    # Calculate bounds for y-level displacement
    bounds_df <- tx_df %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarise(tx_start = min(start), tx_end = max(end)) %>%
      dplyr::arrange(tx_start)

    tx_iranges <- IRanges::IRanges(start = bounds_df$tx_start, end = bounds_df$tx_end)
    bounds_df$y_level <- IRanges::disjointBins(tx_iranges)

    exons_df <- dplyr::left_join(exons_df, bounds_df[, c("SYMBOL", "y_level")], by = "SYMBOL")
    if(nrow(introns_df) > 0) introns_df <- dplyr::left_join(introns_df, bounds_df[, c("SYMBOL", "y_level")], by = "SYMBOL")

    labels_df <- bounds_df %>% dplyr::rename(label_name = SYMBOL)

  } else {
    exons_df <- as.data.frame(exons_list) %>% dplyr::rename(tx_id = group_name)
    introns_df <- as.data.frame(introns_list)
    if(nrow(introns_df) > 0) introns_df <- introns_df %>% dplyr::rename(tx_id = group_name)

    # Calculate bounds for y-level displacement
    bounds_df <- tx_df %>%
      dplyr::rename(tx_start = start, tx_end = end) %>%
      dplyr::arrange(tx_start)

    tx_iranges <- IRanges::IRanges(start = bounds_df$tx_start, end = bounds_df$tx_end)
    bounds_df$y_level <- IRanges::disjointBins(tx_iranges)

    exons_df <- dplyr::left_join(exons_df, bounds_df[, c("TXID", "y_level")], by = c("tx_id" = "TXID"))
    if(nrow(introns_df) > 0) introns_df <- dplyr::left_join(introns_df, bounds_df[, c("TXID", "y_level")], by = c("tx_id" = "TXID"))

    labels_df <- bounds_df %>% dplyr::rename(label_name = SYMBOL)
  }

  return(list(exons = exons_df, introns = introns_df, labels = labels_df))
}


# ---- Plotting ----
get_genes_plot <- function(introns, exons, labels, chrom_start, chrom_end, x_axis = FALSE, color = "black") {

  plot <- ggplot2::ggplot()

  if (nrow(introns) > 0) {
    delta <- (chrom_end - chrom_start) * 0.001
    # creates a minor segment in the middle of each intron indicating the direction of the gene
    introns <- introns %>%
      dplyr::mutate(
        mid = (start + end) / 2,
        arrow_start = ifelse(strand == "-", mid + delta, mid - delta),
        arrow_end   = ifelse(strand == "-", mid - delta, mid + delta)
      )

    plot <- plot +
      ggplot2::geom_segment(data = introns,
                            ggplot2::aes(x = start, xend = end, y = y_level, yend = y_level),
                            color = color, linewidth = 0.5) +

      ggplot2::geom_segment(data = introns,
                            ggplot2::aes(x = arrow_start, xend = arrow_end, y = y_level, yend = y_level),
                            color = color, linewidth = 0.5,
                            arrow = grid::arrow(length = grid::unit(0.05, "inches"), type = "closed"))
  }

  if (nrow(exons) > 0) {
    plot <- plot +
      ggplot2::geom_rect(data = exons,
                         ggplot2::aes(xmin = start, xmax = end,
                                      ymin = y_level - 0.2, ymax = y_level + 0.2),
                         fill = color, color = color)
  }

  if (nrow(labels) > 0) {
    labels <- labels %>% dplyr::mutate(label_x = base::pmax(tx_start, chrom_start))

    plot <- plot +
      ggplot2::geom_text(data = labels,
                         ggplot2::aes(x = label_x, y = y_level, label = label_name),
                         hjust = 0, vjust = -1.2, size = 3.5, fontface = "italic")
  }

  plot <- plot +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::coord_cartesian(xlim = c(chrom_start, chrom_end)) +
    ggplot2::scale_y_continuous(expand = c(0.1, 0.5))

  if (!x_axis) {
    plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                  axis.ticks.x = ggplot2::element_blank(),
                                  axis.line.x = ggplot2::element_blank())
  }

  return(plot)
}


# ---- Exposed functions ----
#' @title Plot genes at specified chromosomal region
#' @description Plot genes and transcripts in a specified genomic region.
#' @param region character: string specifying region to plot (accepts chr:start-end and chr-start-end formats).
#' @param txdb TxDb object: Transcript database, e.g., TxDb.Hsapiens.UCSC.hg38.knownGene.
#' @param orgdb OrgDb object: Organism database for gene symbols, e.g., org.Hs.eg.db.
#' @param collapse bool: If TRUE, merges all transcripts of a gene into a single flattened structure.
#' @param fully_in_view bool: If TRUE, removes any genes that extend beyond the specified region.
#' @param x_axis bool: Whether to draw the x-axis.
#' @param color character: Color for the genes.
#' @returns A ggplot2 plot.
#' @export
genes_track <- function(region, txdb, orgdb, collapse = FALSE, fully_in_view = FALSE, x_axis = TRUE, color = "black") {
  region_split <- split_region_string(region)

  gene_data <- fetch_genes(
    chr = region_split[["chr"]],
    chrom_start = region_split[["start"]],
    chrom_end = region_split[["end"]],
    txdb = txdb,
    orgdb = orgdb,
    collapse = collapse,
    fully_in_view = fully_in_view
  )

  plot <- get_genes_plot(
    introns = gene_data[["introns"]],
    exons = gene_data[["exons"]],
    labels = gene_data[["labels"]],
    chrom_start = region_split[["start"]],
    chrom_end = region_split[["end"]],
    x_axis = x_axis,
    color = color
  )

  return(plot)
}
