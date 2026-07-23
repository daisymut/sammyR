#' Import BigWig Files and Re-bin to Common Intervals
#'
#' Imports one or more BigWig files as RleList and computes the binned-average
#' signal over a shared set of bins. Seqlevels are aligned to `bin_list` before
#' binning so tracks from different sources are comparable. Signal type is
#' irrelevant: single-fraction coverage or MLE ratios are handled identically.
#'
#' @param files Character vector of paths to BigWig files.
#' @param bin_list A GRanges of target bins to average over.
#' @param names Character vector of names for the returned list, one per file.
#' @param genome Optional genome tag (e.g. "hg38"); if NULL (default) not set.
#' @param cores Integer for parallel::mclapply. Default 1 (serial); >1
#'   parallelises across files (ignored on Windows).
#' @return A named list of GRanges, each with a `score` metadata column.
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevels seqlevels<- genome<-
#' @importFrom GenomicRanges binnedAverage
#' @importFrom parallel mclapply
#' @export
import_and_rebin__bw <- function(files, bin_list, names, genome = NULL, cores = 1) {
    bin_names <- GenomeInfoDb::seqlevels(bin_list)
    bws <- parallel::mclapply(files, mc.cores = cores, function(file) {
        bwR <- rtracklayer::import(file, format = "BigWig", as = "RleList")
        missing <- setdiff(bin_names, names(bwR))
        if (length(missing))
            stop("seqlevels di bin_list assenti dal BigWig: ", paste(missing, collapse = ", "))
        bwR <- bwR[bin_names]                                   # riordina, order-safe
        bw  <- GenomicRanges::binnedAverage(bins = bin_list, numvar = bwR, varname = "score")
        if (!is.null(genome)) GenomeInfoDb::genome(bw) <- genome
        bw
  })
  names(bws) <- names
  bws
}

import_and_rebin__bw <- function(files, bin_list, names, genome = NULL, cores = 1) {
  bws <- parallel::mclapply(files, mc.cores = cores, function(file) {
    bwR <- rtracklayer::import(file, format = "BigWig", as = "RleList")

    bin_names <- GenomeInfoDb::seqlevels(bin_list)
    bw_names  <- names(bwR)
    ## exact, order-safe match (replaces the fragile regex grep)
    chr_order <- match(bw_names, bin_names)
    if (anyNA(chr_order))
      stop("BigWig seqlevels not in bin_list: ",
           paste(bw_names[is.na(chr_order)], collapse = ", "))

    bins_for_bw <- bin_list
    GenomeInfoDb::seqlevels(bins_for_bw) <- bin_names[chr_order]

    bw <- GenomicRanges::binnedAverage(
      bins    = bins_for_bw,
      numvar  = bwR[GenomeInfoDb::seqlevels(bins_for_bw)],
      varname = "score"
    )

    if (!is.null(genome)) GenomeInfoDb::genome(bw) <- genome
    bw
  })
  names(bws) <- names
  bws
}

#' Deterministic k-means wrapper (seeded)
#'
#' Thin wrapper around \code{stats::kmeans} with a fixed seed and high
#' \code{nstart}, so subcompartment clustering is reproducible run to run.
#' Adapted from CALDER's internal helper.
#'
#' @param iter.max Maximum iterations passed to \code{stats::kmeans}.
#' @param nstart Number of random starts passed to \code{stats::kmeans}.
#' @param ... Further arguments forwarded to \code{stats::kmeans}
#'   (e.g. \code{x}, \code{centers}).
#'
#' @return A \code{kmeans} object.
#'
#' @importFrom stats kmeans
#' @keywords internal
seeded_kmeans <- function( iter.max = 1000, nstart = 50, ... ){
  set.seed( 1 )
  kmeans( iter.max = iter.max, nstart = nstart, ... )
}

#' Create a TxDb from a GTF Across GenomicFeatures Versions
#'
#' Wrapper around `makeTxDbFromGFF` that works whether the constructor lives in
#' `txdbmaker` (Bioconductor >= 3.19 / GenomicFeatures >= 1.56) or still in
#' `GenomicFeatures` (older pipeline containers). Prefers `txdbmaker` when
#' available and falls back to `GenomicFeatures` otherwise.
#'
#' @param gtf_file Path to a GTF/GFF annotation file.
#' @return A TxDb object.
#' @keywords internal
make_txdb_from_gff <- function(gtf_file) {
  pkg <- if (requireNamespace("txdbmaker", quietly = TRUE)) "txdbmaker" else "GenomicFeatures"
  maker <- getExportedValue(pkg, "makeTxDbFromGFF")
  maker(gtf_file)
}

#' Build a Protein-Coding Gene Annotation from a GTF
#'
#' Constructs a TxDb from a GTF, extracts gene ranges, keeps only genes with at
#' least one protein-coding transcript, and strips version suffixes from gene IDs.
#'
#' @param gtf_file Path to a GTF/GFF annotation file.
#' @return A GRanges of protein-coding genes with cleaned `gene_id` mcols.
#' @importFrom GenomicFeatures genes
#' @importFrom S4Vectors mcols mcols<-
#' @export
setup_gene_annotation <- function(gtf_file) {
  # Create TxDb from GTF
  txdb <- make_txdb_from_gff(gtf_file)
  genes <- genes(txdb)
  # Get protein coding genes
  geneid_codingdf <- summarizeProteinCodingGenes(txdb)
  final_genes <- genes[genes$gene_id %in% geneid_codingdf[geneid_codingdf$n_coding > 0,]$gene]
  # Clean gene IDs
  mcols(final_genes)$gene_id <- gsub("\\..*", "", mcols(final_genes)$gene_id)
  return(final_genes)
}

#' Riassume il contenuto protein-coding per gene
#'
#' Per ogni gene di un TxDb conta i trascritti totali e quanti sono
#' protein-coding (cioè hanno almeno un CDS), restituendo un sommario per gene.
#'
#' @param txdb Un oggetto \code{TxDb}.
#' @return Un data.frame con una riga per gene e colonne \code{gene},
#'   \code{n_tx}, \code{n_coding}, \code{n_non_coding}.
#' @importFrom GenomicFeatures cdsBy transcripts
#' @importFrom S4Vectors mcols splitAsList
#' @importFrom methods is
#' @keywords internal
summarizeProteinCodingGenes <- function(txdb) {
  stopifnot(is(txdb, "TxDb"))
  protein_coding_tx <- names(cdsBy(txdb, use.names = TRUE))
  all_tx <- mcols(transcripts(txdb, columns = c("gene_id", "tx_name")))
  all_tx$gene_id <- as.character(all_tx$gene_id)
  all_tx$is_coding <- all_tx$tx_name %in% protein_coding_tx
  tmp <- splitAsList(all_tx$is_coding, all_tx$gene_id)
  gene <- names(tmp)
  n_tx <- lengths(tmp)
  n_coding <- sum(tmp)
  n_non_coding <- n_tx - n_coding
  data.frame(gene, n_tx, n_coding, n_non_coding, stringsAsFactors = FALSE)
}

#' Add comparison metadata columns to a data frame
#'
#' Tags a data frame with the comparison it belongs to, so rows from many
#' comparisons can be row-bound later and remain traceable.
#'
#' @param df A data frame (typically bins or genes for one comparison).
#' @param comparison_name Character. Name of the comparison.
#' @param current_ratio Character. The ratio/contrast identifier.
#' @param fraction Character. Fraction label (e.g. "all" or a fraction name).
#' @param direction Character. Direction label (e.g. "up", "down", "all").
#'
#' @return The input `df` with `comparison`, `ratio`, `fraction` and
#'   `direction` columns added.
#' @export
add_metadata <- function(df, comparison_name, current_ratio, fraction, direction) {
  df[['comparison']] <- comparison_name
  df[['ratio']] <- current_ratio
  df[['fraction']] <- fraction
  df[['direction']] <- direction
  return(df)
}

#' Save bins data to a structured CSV
#'
#' Row-binds a list of per-category bin data frames, keeps the relevant columns
#' depending on whether this is an "all bins" or a filtered export, and writes a
#' CSV. If the list is empty, nothing is written and a message is printed.
#'
#' @param data_list List of data frames to combine. If empty, nothing is written.
#' @param current_ratio Character. Ratio/contrast identifier, used in the filename.
#' @param comparison_name Character. Comparison name, used in the filename.
#' @param file_suffix Character. Suffix selecting the export type; if it matches
#'   "all_bins", the `fraction`/`direction` columns are dropped.
#' @param g1,g2 Optional character. Reference and test group names; when both are
#'   given, their per-group statistic columns are kept in filtered exports.
#'
#' @return Invisibly `NULL`; called for its side effect of writing a CSV.
#' @importFrom utils write.csv
#' @export
save_bins_data <- function(data_list, current_ratio, comparison_name, file_suffix, g1 = NULL, g2 = NULL) {
  if (length(data_list)) {
    df <- do.call(rbind, data_list)
    if (grepl("all_bins", file_suffix, ignore.case = TRUE)) {
      cols_to_remove <- c('fraction', 'direction')
      available_cols <- colnames(df)
      cols_to_keep <- setdiff(available_cols, cols_to_remove)
      df <- df[, cols_to_keep]
    } else {
      base_cols <- c('seqnames', 'start', 'end', 'ratio', 'comparison', 'fraction', 'direction')
      essential_stats_cols <- c()
      if (!is.null(g1) && !is.null(g2)) {
        essential_stats_cols <- c(
          paste0(g1, "_serrx2_lower"), paste0(g1, "_serrx2_upper"),
          paste0(g1, "_mean"), paste0(g1, "_serrX2"),
          paste0(g2, "_serrx2_lower"), paste0(g2, "_serrx2_upper"),
          paste0(g2, "_mean"), paste0(g2, "_serrX2"),
          "delta",
          "cohen.estimate",
          "cohen.magnitude"
        )}
      essential_cols <- c(base_cols, essential_stats_cols)
      available_cols <- colnames(df)
      cols_to_keep <- intersect(essential_cols, available_cols)
      df <- df[, cols_to_keep]
    }
    output_file <- paste0(current_ratio, "_", comparison_name, "_", file_suffix, ".csv")
    write.csv(df, file = output_file, quote = FALSE, row.names = FALSE)
  } else {
    cat("No", file_suffix, "data found for", comparison_name, "\\n")
  }
  invisible( NULL )
}
