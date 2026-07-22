#' Select and Classify Differentially Soluble Bins
#'
#' Core engine of the differential solubility analysis. For a single group
#' comparison, computes per-bin confidence intervals for each group, the delta
#' between them, directional shift calls and Cohen's d, then classifies bins
#' into up/down categories for two fractions and, when annotation is available,
#' maps them onto overlapping gene promoters.
#'
#' @param combination Character scalar naming the comparison; used as the column
#'   selector into the global `pr` matrix to retrieve the two group labels.
#' @param allmixeddf_grobj A GRanges of bins with per-sample score columns in
#'   `mcols` (typically quantile-normalised).
#' @param fraction1 Character label of the first fraction (default "S2S").
#' @param fraction2 Character label of the second fraction (default "S3").
#' @param ths Numeric absolute-mean threshold for retaining a bin (default 0.1).
#'
#' @return A named list keyed by group labels, containing:
#'   \itemize{
#'     \item \code{<test>_vs_<ref>_all_shifting_bins} - all classified shifting bins;
#'     \item four per-fraction category GRanges (\code{<frac>_up} / \code{<frac>_down});
#'     \item \code{<test>_genes_<ref>} - gene counts per category;
#'     \item \code{genes} - vectors of overlapping gene IDs per category;
#'     \item \code{<test>_gr_<ref>} - all threshold-passing bins;
#'     \item \code{<test>_allgr_<ref>} - all bins with Cohen's d.
#'   }
#'
#' @note Currently relies on `pr` (group matrix) and, for gene mapping,
#'   `final_genes` being present in `.GlobalEnv`. This global-state dependency is
#'   slated for refactor into explicit arguments.
#'
#' @seealso \code{\link{import_and_rebin__bw}}
#'
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps promoters
#' @export
Bins_selector <- function(combination, allmixeddf_grobj, fraction1 = "S2S", fraction2 = "S3", ths = 0.1) {
  cat("Running Bins_selector for combination:", combination, "\n")

  # Get pr from global environment
  pr <- get("pr", envir = .GlobalEnv)

  # Define groups to compare and select their samples
  x <- get(pr[, combination][1])
  y <- get(pr[, combination][2])
  xgroup <- pr[, combination][1]
  ygroup <- pr[, combination][2]

  cat("Groups:", paste(x, collapse = ", "), "vs", paste(y, collapse = ", "), "\n")

  new_selection <- allmixeddf_grobj
  mcols(new_selection) <- mcols(new_selection)[c(x, y)]

  # Calculations
  confint_mean_serr_first <- apply(as.matrix(mcols(new_selection)[c(x)]), 1, confidence_interval, nm = xgroup)
  confint_mean_serr_second <- apply(as.matrix(mcols(new_selection)[c(y)]), 1, confidence_interval, nm = ygroup)
  delta <- confint_mean_serr_first[1, ] - confint_mean_serr_second[1, ]
  res <- cbind(t(confint_mean_serr_first), t(confint_mean_serr_second), delta)
  df_toadd1 <- do.call("cbind", as.data.frame(res))
  mcols(new_selection) <- cbind(mcols(new_selection), df_toadd1)

  # Range analysis forward comparison
  range_analysis <- lapply(1:length(y), function(n) {
    y_name <- y[n]
    z <- apply(as.matrix(mcols(new_selection)[c(paste0(xgroup, "_serrx2_lower"), paste0(xgroup, "_serrx2_upper"), paste0(xgroup, "_mean"), y_name)]), 1,
               is_in_serrx2_range_and_shift,
               xgroup_name = xgroup,
               y_name = y_name
    )
    return(as.data.frame(t(z)))
  })

  df_toadd <- do.call("cbind", range_analysis)
  mcols(new_selection) <- cbind(mcols(new_selection), df_toadd)

  # Range analysis reverse comparison
  range_analysis_rev <- lapply(1:length(x), function(n) {
    x_name <- x[n]
    z <- apply(as.matrix(mcols(new_selection)[c(paste0(ygroup, "_serrx2_lower"), paste0(ygroup, "_serrx2_upper"), paste0(ygroup, "_mean"), x_name)]), 1,
               is_in_serrx2_range_and_shift,
               xgroup_name = ygroup,
               y_name = x_name
    )
    return(as.data.frame(t(z)))
  })

  df_toadd2 <- do.call("cbind", range_analysis_rev)
  mcols(new_selection) <- cbind(mcols(new_selection), df_toadd2)

  # Select bins
  prvdf <- as.data.frame(new_selection)
  prvdf_gr <- makeGRangesFromDataFrame(prvdf, keep.extra.columns = TRUE)

  # Calculate Cohen's d for all bins
  all_bins_with_cohens <- func_ztest_gr_byrow(prvdf_gr,
                                              x = x,
                                              y = y,
                                              cohenthresh = 0  # No filtering by Cohen's d here
  )
  cat("Cohen's d calculated for", length(all_bins_with_cohens), "bins\n")

  prvdf_with_cohens <- as.data.frame(all_bins_with_cohens)
  prvdftest <- prvdf_with_cohens[abs(prvdf_with_cohens[paste0(xgroup, "_mean")]) >= ths, ]
  cat("After threshold filtering:", nrow(prvdftest), "bins remain\n")

  pprvlow <- prvdftest[paste0(y, "_shift")] == "lower"   # <-- cambiato da _ov_specs a _shift
  pprvhigh <- prvdftest[paste0(y, "_shift")] == "higher" # <-- cambiato da _ov_specs a _shift

  # Sum up by group 1
  prvdftest[, paste0(ygroup, "_sign_SUM")] <- rowSums(sapply(prvdftest[, paste0(y, "_sign")], as.numeric)) # cambiato da _ov_check a _sign
  prvdftest$ovlow <- apply(pprvlow, 1, sum) * -1
  prvdftest$ovvhigh <- apply(pprvhigh, 1, sum)

  # Commutative group testing
  pprvlow_X <- prvdftest[paste0(x, "_shift")] == "lower"   #  <-- cambiato da  ov_specs a _shift
  pprvhigh_X <- prvdftest[paste0(x, "_shift")] == "higher" #  <-- cambiato da  ov_specs a _shift

  # Sum up by group 2
  prvdftest[, paste0(xgroup, "_sign_SUM")] <- rowSums(sapply(prvdftest[, paste0(x, "_sign")], as.numeric)) ## <-- cambiato da _ov_check a _sign

  # Count characteristics
  prvdftest$ovlow_X <- apply(pprvlow_X, 1, sum) * -1
  prvdftest$ovvhigh_X <- apply(pprvhigh_X, 1, sum)

  # Add meantosep column
  prvdftest$meantosep <- prvdftest[,paste0(xgroup, "_mean_startsign")]

  prvdftest_gr <- makeGRangesFromDataFrame(prvdftest, keep.extra.columns = TRUE)

  # Apply statistical testing
  up_down_to_ztest_grr <- prvdftest_gr[abs(mcols(prvdftest_gr)$cohen.estimate) >= 2] # huge

  cat("Statistical testing completed. Regions passing threshold:", length(up_down_to_ztest_grr), "\n")

  # Separate bins according to meantosep
  startmeanpos <- up_down_to_ztest_grr[up_down_to_ztest_grr$meantosep >= 0]
  startmeanneg <- up_down_to_ztest_grr[up_down_to_ztest_grr$meantosep <= 0]

  # Select coherent bins based on direction of change
  ovvhighconservedpos <- startmeanpos[abs(startmeanpos$ovvhigh) > abs(startmeanpos$ovlow)]
  ovlowconservedpos <- startmeanpos[abs(startmeanpos$ovvhigh) < abs(startmeanpos$ovlow)]
  ovvhighconservedneg <- startmeanneg[abs(startmeanneg$ovvhigh) > abs(startmeanneg$ovlow)]
  ovlowconservedneg <- startmeanneg[abs(startmeanneg$ovvhigh) < abs(startmeanneg$ovlow)]

  # Add bin type labels
  ovvhighconservedpos$bintype <- rep(paste0(fraction1, "_up"), length(ovvhighconservedpos))
  ovlowconservedpos$bintype <- rep(paste0(fraction1, "_down"), length(ovlowconservedpos))
  ovvhighconservedneg$bintype <- rep(paste0(fraction2, "_up"), length(ovvhighconservedneg))
  ovlowconservedneg$bintype <- rep(paste0(fraction2, "_down"), length(ovlowconservedneg))

  # Control groups
  controlstartmeanpos <- prvdftest_gr[prvdftest_gr$meantosep >= 0]
  controlstartmeanneg <- prvdftest_gr[prvdftest_gr$meantosep <= 0]

  controlstartmeanpos$bintype <- rep(fraction1, length(controlstartmeanpos))
  controlstartmeanneg$bintype <- rep(fraction2, length(controlstartmeanneg))

  controlstartmeanpos <- controlstartmeanpos[!controlstartmeanpos %in%
                                               c(ovvhighconservedpos, ovlowconservedpos)]
  controlstartmeanneg <- controlstartmeanneg[!controlstartmeanneg %in%
                                               c(ovvhighconservedneg, ovlowconservedneg)]

  # All groups to return
  all_gr_toreturn <- c(ovvhighconservedpos,
                       ovlowconservedpos,
                       ovvhighconservedneg,
                       ovlowconservedneg,
                       controlstartmeanpos,
                       controlstartmeanneg)

  # Make a list of bins to save and analyze
  list_ofbins_to_save_and_analyse <- setNames(
    list(
      ovvhighconservedpos,
      ovlowconservedpos,
      ovvhighconservedneg,
      ovlowconservedneg
    ),
    c(
      paste0(fraction1, "_up"),
      paste0(fraction1, "_down"),
      paste0(fraction2, "_up"),
      paste0(fraction2, "_down")
    )
  )

  # Gene analysis (if gtf provided)
  list_of_vector_geneNumber <- list()
  list_of_genes_vec <- list()

  final_genes_obj <- tryCatch(get("final_genes", envir = .GlobalEnv), error = function(e) NULL)

  if (!is.null(final_genes_obj)) {
    for (i in names(list_ofbins_to_save_and_analyse)) {
      gr_touse <- list_ofbins_to_save_and_analyse[[i]]
      if (length(gr_touse) != 0) {
        cat(paste0(ygroup, "_vs_", xgroup, "_", i, " has ", length(gr_touse), " regions"), "\n")

        # Calculate overlapping genes
        genes_gr <- final_genes_obj[findOverlaps(gr_touse, promoters(final_genes_obj, upstream = 2500, downstream = 500))@to]

        list_of_vector_geneNumber[[i]] <- length(mcols(genes_gr)$gene_id)
        list_of_genes_vec[[paste0(ygroup, "_vs_", xgroup, "_", i)]] <- unique(mcols(genes_gr)$gene_id)
      } else {
        list_of_vector_geneNumber[[i]] <- 0
        list_of_genes_vec[[paste0(ygroup, "_vs_", xgroup, "_", i)]] <- character(0)
      }
    }
  } else {
    for (i in names(list_ofbins_to_save_and_analyse)) {
      list_of_vector_geneNumber[[i]] <- 0
      list_of_genes_vec[[paste0(ygroup, "_vs_", xgroup, "_", i)]] <- character(0)
    }
  }

  # Numeric coding for the groups
  mcols(ovvhighconservedpos)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(2, length(ovvhighconservedpos))
  mcols(ovlowconservedpos)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(1, length(ovlowconservedpos))
  mcols(ovvhighconservedneg)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(-1, length(ovvhighconservedneg))
  mcols(ovlowconservedneg)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(-2, length(ovlowconservedneg))

  # Prepare results list
  x <- list()
  x[[paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins")]] <- c(
    ovvhighconservedpos,
    ovlowconservedpos,
    ovvhighconservedneg,
    ovlowconservedneg
  )
  x[[names(list_ofbins_to_save_and_analyse[1])]] <- list_ofbins_to_save_and_analyse[[1]]
  x[[names(list_ofbins_to_save_and_analyse[2])]] <- list_ofbins_to_save_and_analyse[[2]]
  x[[names(list_ofbins_to_save_and_analyse[3])]] <- list_ofbins_to_save_and_analyse[[3]]
  x[[names(list_ofbins_to_save_and_analyse[4])]] <- list_ofbins_to_save_and_analyse[[4]]
  x[[paste0(ygroup, "_genes_", xgroup)]] <- list_of_vector_geneNumber
  x[["genes"]] <- list_of_genes_vec
  x[[paste0(ygroup, "_gr_", xgroup)]] <- prvdftest_gr
  x[[paste0(ygroup, "_allgr_", xgroup)]] <- all_bins_with_cohens

  names(x) <- c(
    paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins"),
    paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[1]), "_", xgroup),
    paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[2]), "_", xgroup),
    paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[3]), "_", xgroup),
    paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[4]), "_", xgroup),
    paste0(ygroup, "_genes_", xgroup),
    "genes",
    paste0(ygroup, "_gr_", xgroup),
    paste0(ygroup, "_allgr_", xgroup)
  )

  cat("Results:", paste0(names(x[paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins")]), "_", length(x[[paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins")]]), "_bins"), "\n")
  return(x)
}

#' Confidence Interval via Twice the Standard Error
#'
#' Computes a mean-centred interval of +/- two standard errors for a numeric
#' vector, returning bounds, mean and 2xSE width with prefixed names.
#'
#' @param vector Numeric vector of observations.
#' @param nm Character prefix for the returned names (default "prove").
#' @return Named numeric vector: `<nm>_serrx2_lower`, `<nm>_serrx2_upper`,
#'   `<nm>_mean`, `<nm>_serrX2`.
#' @importFrom stats sd
#' @keywords internal
confidence_interval <- function(vector, nm="prove") {
  # Standard deviation of sample
  vec_serr <- sd(vector)/sqrt(length(vector))
  vec_serr2x<- vec_serr*2
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  nm_confint_low<-paste0(nm,"_serrx2_lower")
  nm_confint_up<-paste0(nm,"_serrx2_upper")
  name_mean <- paste0(nm,"_mean")
  name_serr <- paste0(nm,"_serrX2")
  result <- c(nm_confint_low = vec_mean - vec_serr2x,
              nm_confint_up = vec_mean + vec_serr2x,
              name_mean = vec_mean,
              name_serr = vec_serr2x
  )
  names(result) <- c(nm_confint_low,nm_confint_up,name_mean,name_serr)
  return(result)
}

#' Locate a Group Value Relative to Another Group's 2xSE Interval
#'
#' Tests whether a group's value falls inside, below or above the twice-standard-
#' error interval of a reference group, and derives the sign, direction and
#' solubility-shift labels for that bin.
#'
#' @param vector Named numeric vector with the reference group's `_serrx2_lower`,
#'   `_serrx2_upper`, `_mean` entries plus the tested value.
#' @param xgroup_name Character prefix identifying the reference group entries.
#' @param y_name Character name of the tested value entry.
#' @return Named character vector: `<y_name>_sign`, `<y_name>_shift`,
#'   `<y_name>_sol_shift`, `<xgroup_name>_mean_startsign`.
#' @importFrom data.table between
#' @keywords internal
is_in_serrx2_range_and_shift <- function(vector,
                                         xgroup_name,
                                         y_name
){
  #meann,lower,upper,mean_tp_name,whereis_tp_name){
  xgroup_lower_bound_confint_name<-paste0(xgroup_name,"_serrx2_lower")
  xgroup_upper_bound_confint_name<-paste0(xgroup_name,"_serrx2_upper")
  xgroup_mean_name<-paste0(xgroup_name,"_mean")

  xgroup_lower_bound_confint<-vector[xgroup_lower_bound_confint_name][1]
  xgroup_upper_bound_confint<-vector[xgroup_upper_bound_confint_name][1]

  yname_val<-vector[y_name][1]
  xgroup_mean<-vector[xgroup_mean_name][1]

  confint_check <-''
  mean_startsign<- ''
  mean_sign <-''
  shift_solubility<- ''

  col_is_in_confint_name<-paste0(y_name,"_sign") # cambiato da ov_check a _sign
  col_whereis_tp_name<- paste0(y_name,"_shift") # cambiato da ov_specs a _shift
  shift_solubility_name<- paste0(y_name, "_sol_shift")
  mean_startsign_name<- paste0(xgroup_name, "_mean_startsign")
  ##########check if both or one mean is in the confint of the other mean
  if ( between(yname_val, xgroup_lower_bound_confint, xgroup_upper_bound_confint)
  ) {
    confint_check <- sign(yname_val)
    mean_sign <- "nodiff"
    shift_solubility <- "constant_solubility"
    mean_startsign<- sign(xgroup_mean)
    ##########check if confint are not overlapping, the ygroup timepoint confint is lower than xgroup
  }else if (yname_val < xgroup_lower_bound_confint ) {

    confint_check <- sign(yname_val)

    #mean_sign <- paste0(y_name ,"_lower_than_",xgroup_name)
    mean_sign <- "lower"
    shift_solubility <- check_sign(yname_val,xgroup_mean)
    mean_startsign<- sign(xgroup_mean)

    ##########check if confint are not overlapping, the xgroup timepoint confint is lower than ygroup
  }else if (xgroup_upper_bound_confint < yname_val) {

    confint_check <-sign(yname_val)
    mean_sign <- "higher"
    shift_solubility <- check_sign(yname_val,xgroup_mean)
    mean_startsign<- sign(xgroup_mean)
  }

  else{

    confint_check <- sign(yname_val)
    mean_sign <- "no_idea"
    shift_solubility <- check_sign(yname_val,xgroup_mean)
    mean_startsign<- sign(xgroup_mean)
  }

  result <- c(confint_check,mean_sign,shift_solubility,mean_startsign)

  names(result) <- c(col_is_in_confint_name,col_whereis_tp_name,shift_solubility_name,mean_startsign_name)
  return(result)
}

#' Per-Bin Effect Size (Cohen's d) with Magnitude Filtering
#'
#' For each bin (row) in a GRanges, computes Cohen's d between two groups of
#' score columns, appends the estimate and its qualitative magnitude to the
#' metadata, and keeps only bins whose absolute effect size meets the threshold.
#'
#' @param gr A GRanges whose mcols hold the per-sample score columns.
#' @param x Selector for the first group's columns.
#' @param y Selector for the second group's columns.
#' @param cohenthresh Numeric absolute Cohen's d cut-off (default 0.8).
#' @return The input GRanges with `cohen.estimate` and `cohen.magnitude` mcols
#'   added, filtered to |cohen.estimate| >= cohenthresh.
#' @note A BSDA z-test and BH correction are scaffolded but commented out; the
#'   current criterion is effect size only. Depends on the internal `extract_nth`.
#' @importFrom effsize cohen.d
#' @importFrom S4Vectors mcols mcols<-
#' @keywords internal


func_ztest_gr_byrow <- function(gr,
                                x,
                                y,
                                # correction_method = "BH",
                                cohenthresh = 0.8)
{
  ppval <- lapply(seq(nrow(as.data.frame(mcols(gr)))), function(i) {

    # Cohen's d calculation
    cohend <- cohen.d(
      unlist(as.vector(as.data.frame(mcols(gr))[x][i,])),
      unlist(as.vector(as.data.frame(mcols(gr))[y][i,]))
    )
    cohen.estimate <- cohend$estimate
    cohend.magnitude <- as.character(cohend$magnitude)

    # BSDA Z-test implementation
    # ztest <- z.test(x = as.data.frame(mcols(gr))[x][i,],
    #                 y = as.data.frame(mcols(gr))[y][i,],
    #                 sigma.x = sd(as.data.frame(mcols(gr))[x][i,]),
    #                 sigma.y = sd(as.data.frame(mcols(gr))[y][i,]),
    #                 alternative = 'two.sided',
    #                 conf.level = 0.99)

    # ztest_pvalue <- ztest$p.value

    zzzz <- list(cohen.estimate, cohend.magnitude)
    names(zzzz) <- c("cohen.estimate", "cohen.magnitude")
    return(zzzz)
  })

  df_tomerge_mcols <- data.frame(
    cohen.estimate = unlist(extract_nth(ppval, 1)),
    cohen.magnitude = unlist(extract_nth(ppval, 2))
  )

  # Apply Benjamini-Hochberg correction
  # df_tomerge_mcols[[paste0("ztest_", correction_method, "_correct")]] <- p.adjust(df_tomerge_mcols$ztest, method = correction_method)

  mcols(gr) <- cbind(mcols(gr), df_tomerge_mcols)

  # Filter by Cohen's d
  gr <- gr[abs(mcols(gr)$cohen.estimate) >= cohenthresh]
  return(gr)
}

#' Classify Solubility Shift Direction from Sign Change
#'
#' Compares the sign of a value against the sign of a reference mean to label
#' the direction of a solubility shift.
#'
#' @param x Numeric value (e.g. a group mean).
#' @param meann Numeric reference mean to compare against.
#' @return "constant_solubility", "shift_increase" or "shift_decrease".
#' @keywords internal
check_sign <- function(x,meann){
  if ( sign(x) == sign(meann) ){
    return("constant_solubility")
  } else if (sign(x) > sign(meann)) {
    return("shift_increase")
  }else if (sign(x) < sign(meann)) {
    return("shift_decrease")
  }
}

#' Extract the n-th Element from Each Member of a List
#'
#' Convenience extractor equivalent to `sapply(lst, `[`, n)`.
#'
#' @param lst A list whose elements are each indexable by `n`.
#' @param n Index (name or position) to extract from each element.
#' @return A vector of extracted elements, one per list member.
#' @keywords internal
extract_nth <- function(lst, n){
  sapply(lst, `[`, n)
}
