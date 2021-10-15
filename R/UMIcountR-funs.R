# extract_spike_dat -------------------------------------------------------

#' @title Load and extract spUMIs from sequencing reads
#' @description extract_spike_dat is used to extract and parse molecular spike reads from bam files.
#' @param bam_path path to input bam file, must be indexed zUMIs output file
#' @param spikename name of the geneID for molecular spikes, Default: 'g_diySpike4'
#' @param spikecontig name of the reference contig for molecular spikes, Default: 'diySpike'
#' @param spikeUMI_start optional: fixed start position of the spUMI within reads, Default: NULL
#' @param spikeUMI_end optional: fixed end position of the spUMI within reads, Default: NULL
#' @param fixed_start_pos optional: fixed start position of the spUMI within reads, Default: NULL
#' @param match_seq_before_UMI sequence to match before the spUMI, Default: NULL
#' @param match_seq_after_UMI sequence to match after the spUMI, Default: NULL
#' @return returns a data.table with reads, their UMI and the raw & error-corrected spUMI.
#' @details The spUMI can be extracted by known position within reads mapping to the molecular spike, when reads are required to map at a specific location (eg. for Smart-seq3 5' reads). In this case, all three arguments spikeUMI_start, spikeUMI_end and fixed_start_pos are required. When the position of the spUMI within reads is variable, provide the parameters match_seq_before_UMI and match_seq_after_UMI to extract the spUMI by matching the known surrounding sequence.
#' @examples 
#' \dontrun{
#' example_dat <- extract_spike_dat(
#'  bam_path = bam,
#'  match_seq_before_UMI = "GAGCCTGGGGGAACAGGTAGG",
#'  match_seq_after_UMI = "CTCGGAGGAGAAA",
#'  spikename = "g_diySpike4", spikecontig = "diySpike",
#' )
#' }
#' @seealso 
#'  \code{\link[Rsamtools]{BamInput}},\code{\link[Rsamtools]{ScanBamParam}}
#'  \code{\link[data.table]{data.table-package}}
#' @rdname extract_spike_dat
#' @export 
#' @importFrom Rsamtools idxstatsBam ScanBamParam scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @import data.table
#' 
extract_spike_dat <- function(bam_path, spikename = "g_diySpike4", spikecontig = "diySpike",
                              spikeUMI_start = NULL , spikeUMI_end = NULL, fixed_start_pos = NULL,
                              match_seq_before_UMI = NULL, match_seq_after_UMI = NULL, spikeUMI_length = NULL){
  idxst <- Rsamtools::idxstatsBam(bam_path)
  reflen <- idxst[idxst$seqnames==spikecontig,]$seqlength
  taglist <- c("BC", "QU", "UX", "UB","GE")
  whatlist <- c("rname","pos","cigar","seq")
  
  parms <- Rsamtools::ScanBamParam(tag=taglist,
                                   what=whatlist,
                                   tagFilter = list(GE = spikename),
                                   which = GenomicRanges::GRanges(seqnames = spikecontig, ranges = IRanges::IRanges(1,reflen)))
  print("Reading in data from bam file...")
  dat <- Rsamtools::scanBam(file = bam_path, param = parms)
  dat <- data.table::data.table(
    contig = dat[[1]]$rname,
    pos = dat[[1]]$pos,
    CIGAR = dat[[1]]$cigar,
    seq = as.character(dat[[1]]$seq),
    BC = dat[[1]]$tag$BC,
    QU = dat[[1]]$tag$QU,
    UX = dat[[1]]$tag$UX,
    UB = dat[[1]]$tag$UB
  )
  dat <- dat[!UX==""] #escape internal reads for Smart-seq3
  
  if(!is.null(fixed_start_pos)){
    dat <- dat[pos == fixed_start_pos]
  }
  
  if(is.null(spikeUMI_start) & is.null(match_seq_before_UMI)){
    stop("give either spikeUMI_start or match_seq_before_UMI")
  }
  
  
  if(!is.null(match_seq_before_UMI) & !is.null(match_seq_after_UMI)){
    
    dat[, c("TSSseq","tmpseq") := tstrsplit(seq, match_seq_before_UMI, keep = 1:2)][
      !is.na(tmpseq), c("spikeUMI","seqAfterUMI") := tstrsplit(tmpseq, match_seq_after_UMI, keep = 1:2)][
        , tmpseq := NULL][
          , TSSseq := paste0(TSSseq,match_seq_before_UMI)][
            , seqAfterUMI := paste0(match_seq_after_UMI,seqAfterUMI)]
    
  }else{
    if(is.null(fixed_start_pos) & !is.null(spikeUMI_start)){
      print("Warning: using spike extraction by position in read without fixed read start position!")
    }
    dat[,TSSseq := substr(seq,1,spikeUMI_start-1)][
        ,spikeUMI := substr(seq,spikeUMI_start,spikeUMI_end) ][
        ,seqAfterUMI := substr(seq,spikeUMI_end+1,48)]
  }
  
  dat <- dat[!is.na(spikeUMI)]
  if(!is.null(spikeUMI_length)){
    dat <- dat[nchar(spikeUMI) == spikeUMI_length]
  }
  
  print("Hamming correct spikeUMIs...")
  dat[, spikeUMI_hd1 := UMIcountR::return_corrected_umi(spikeUMI, editham = 1), by = "BC"][
      , spikeUMI_hd2 := UMIcountR::return_corrected_umi(spikeUMI, editham = 2), by = "BC"]
  
  return(dat)
}


# return_corrected_umi ----------------------------------------------------

#' @title Return corrected UMI sequences
#' @description return_corrected_umi performs error correction on an input vector of raw UMI sequences.
#' @param umi_input input character vector of uncorrected UMI sequences
#' @param editham edit distance (hamming) used for collapse, Default: 1
#' @param collapse_mode collapse mode to use, Default: adjacency
#' @return returns a vector of error-corrected UMIs with same length as input sequences. Implemented collapse algorithms: "adjacency","adjacency_directional","adjacency_singleton","cluster"
#' @details For a description of implemented collapse algorithms, please refer to Ziegenhain, Hendriks et al., 2021.
#' @examples 
#' \dontrun{
#' return_corrected_umi(UX_strings, editham = 1, collapse_mode = "adjacency")
#' }
#' @seealso 
#'  \code{\link[data.table]{data.table-package}}
#' @rdname return_corrected_umi
#' @export 
#' @import data.table
#' @import reticulate
#'
return_corrected_umi <- function(umi_input, editham = 1, collapse_mode = NULL){
  if(is.null(collapse_mode)) collapse_mode = "adjacency"
  if(!collapse_mode %in% c("adjacency","adjacency_directional","adjacency_singleton","cluster")) stop("incorrect collapse_mode")
  
  ham_maps <- .hammingFilter_bktree(umiseq = umi_input, edit = editham, collapse_mode = collapse_mode)
  
  if(!"falseUMI" %in% colnames(ham_maps)){
    return(umi_input)
  }else{
    setkey(ham_maps, falseUMI)
    out_dt <- data.table(raw = umi_input)
    out_dt[,ham := raw]
    out_dt[raw %in% ham_maps$falseUMI, ham := ham_maps[out_dt[raw %in% ham_maps$falseUMI]$raw]$trueUMI]
    return(out_dt$ham)
  }
}



# subsample_recompute -----------------------------------------------------
#' @title Subsample molecular spikes per cell and recompute hamming distance correction
#' @description return_corrected_umi performs error correction on an input vector of raw UMI sequences.
#' @param dat input data.table of loaded spUMI data. Must have BC,UX and spikeUMI_hd2 columns.
#' @param mu_nSpikeUMI mean expression level of the spUMI to subsample to.
#' @param threads number of threads to use, Default: 8
#' @return returns ...
#' @details For each cell barcode, subsample the number of spUMIs to an expression level of the given mean and recomputes the error corrected UMIs.
#' @examples 
#' \dontrun{
#' subsample_recompute(spikedat, mu_nSpikeUMI)
#' }
#' @seealso 
#'  \code{\link[data.table]{data.table-package}}
#' @rdname subsample_recompute
#' @export 
#' @import data.table
#' @importFrom parallel mclapply
#'
subsample_recompute <- function(dat, mu_nSpikeUMI, threads = 8 ){
  setDTthreads(1, restore_after_fork = FALSE)
  nspikes_per_bc <- dat[,.(nspikes = length(unique(spikeUMI_hd2))), by = BC]
  dat <- dat[BC %in% nspikes_per_bc[nspikes >= 0.8*mu_nSpikeUMI]$BC]
  dat_l <- split(dat, by = "BC")
  out_l <- parallel::mclapply(dat_l, function(x){
    nspikes <- abs( round( rnorm(mean=mu_nSpikeUMI,sd=sqrt(mu_nSpikeUMI),n=1) ) )
    if(nspikes < length(unique(x$spikeUMI_hd2))){
      spikeIDs <- sample( x = unique(x$spikeUMI_hd2),size = nspikes, replace = FALSE )
      x_out <- x[spikeUMI_hd2 %in% spikeIDs, c("BC","spikeUMI_hd2","UX") , with = FALSE]
    }else{
      x_out <- x[, c("BC","spikeUMI_hd2","UX") , with = FALSE]
    }
    x_out[,UB_hd1 := return_corrected_umi(UX, editham = 1), by ="BC"][
          ,UB_hd2 := return_corrected_umi(UX, editham = 2), by ="BC"]
    return(x_out)
  }, mc.cores = threads)
  out <- rbindlist(out_l)
  out[, mean_nSpikeUMI := mu_nSpikeUMI]
  return(out)
}



# get_overrepresented_spikes ----------------------------------------------

#' @title Find overrepresented spUMI sequences
#' @description get_overrepresented_spikes is used to find the spUMIs of overrepresented molecular spikes.
#' @param dat input data.table of loaded spUMI data. Must have BC and spikeUMI_hd2 columns.
#' @param readcutoff maximum number of raw reads a spUMI should be observed in, Default: 100
#' @param nbccutoff maximum number of cell barcodes a spUMI should be observed in, Default: 5
#' @return returns a list of length 3, where the first list element "plots" contains descriptive plots of chosen filtering cutoffs and the list elements "over_readcutoff" and "over_nbcs" contain the spUMI identities that should be discarded.
#' @details Visualise filtering plots by accessing the "plots" list element of the output. over_readcutoff = spUMIs that are more abundant than the set cutoff of number of reads. over_nbcs = spUMIs that are observed in more cell barcodes than the set cutoff.
#' @examples 
#' \dontrun{
#' get_overrepresented_spikes(spikedat, readcutoff = 100, nbccutoff = 5)
#' }
#' @seealso 
#'  \code{\link[data.table]{data.table-package}}
#' @rdname get_overrepresented_spikes
#' @export 
#' @import data.table
#' @importFrom ggplot2 ggplot theme_classic xlab ylab geom_vline geom_point
#'
get_overrepresented_spikes <- function(dat, readcutoff = 100, nbccutoff = 5){
  userastr <- FALSE
  if("ggrastr" %in% rownames(installed.packages())){
    requireNamespace("ggrastr", quietly = TRUE)
    userastr <- TRUE
  }
  spikeoccurance <- dat[,.(.N, nBCs = length(unique(BC))), by = spikeUMI_hd2]
  setorder(spikeoccurance, N)
  spikeoccurance[,cs := cumsum(N)]
  
  p1 <- ggplot2::ggplot(spikeoccurance,ggplot2::aes(x = N, y = cs)) + 
    ggplot2::theme_classic() + 
    ggplot2::xlab("Reads per spike UMI") + 
    ggplot2::ylab("Cumulative Reads") + 
    ggplot2::geom_vline(xintercept = readcutoff, linetype = "dashed")
  
  setorder(spikeoccurance, nBCs)
  spikeoccurance[,cs2 := (seq(.N))]
  p2 <- ggplot2::ggplot(spikeoccurance,ggplot2::aes(x = nBCs, y = cs2)) + 
    ggplot2::theme_classic() + 
    ggplot2::xlab("Spike in n BCs") + 
    ggplot2::ylab("Cumulative Spikes") + 
    ggplot2::geom_vline(xintercept = nbccutoff, linetype = "dashed")
  
  if(userastr){
    p1 <- p1 + ggrastr::geom_point_rast()
    p2 <- p2 + ggrastr::geom_point_rast()
  }else{
    p1 <- p1 + ggplot2::geom_point()
    p2 <- p2 + ggplot2::geom_point()
  }
  
  overrep1 <- spikeoccurance[N>readcutoff]$spikeUMI
  overrep2 <- spikeoccurance[nBCs > nbccutoff]$spikeUMI
  return(list(plots = list(p1,p2),over_readcutoff = overrep1, over_nbcs = overrep2))
}


# plot_spike_distances ----------------------------------------------------
#' @title Plot minimal pairwise distances between spUMIs within / over cells
#' @description Plot minimal pairwise distances between spUMIs within and over cells to determine the appropriate error-correction edit distance for spUMI sequences
#' @param dat input data.table of loaded spUMI data. Must have BC and spikeUMI columns.
#' @param threads number of threads to use for edit distance calculation, Default: 32
#' @return returns a ggplot object
#' @details uses hamming distance as the edit distance measure.
#' @examples 
#' \dontrun{
#' plot_spike_distances(spikedat)
#' }
#' @seealso 
#'  \code{\link[data.table]{data.table-package}}
#' @rdname plot_spike_distances
#' @export 
#' @import data.table
#' @importFrom stringdist stringdistmatrix
#' @importFrom matrixStats rowMins
#' @importFrom ggplot2 ggplot theme_classic xlab ylab geom_vline geom_point
#'
plot_spike_distances <- function(dat, threads = 32){
  spike_dist_dat <- unique(dat[, c("spikeUMI","BC"), with = F])
  
  spike_dist_dat_out <- lapply(unique(spike_dist_dat$BC), function(b){
    thisBC_spikes <- spike_dist_dat[BC == b]$spikeUMI
    within_dat <- stringdist::stringdistmatrix(thisBC_spikes,thisBC_spikes, method = 'hamming', nthread = threads)
    within_dat[upper.tri(within_dat, diag = TRUE)] <- NA
    within_mins <- matrixStats::rowMins(within_dat, na.rm = T)
    between_dat <- stringdist::stringdistmatrix(thisBC_spikes, 
                                                spike_dist_dat[!BC == b][,.SD[sample(.N,length(thisBC_spikes), replace = F)]]$spikeUMI,
                                                method= 'hamming', nthread = threads)
    between_dat[upper.tri(between_dat, diag = TRUE)] <- NA
    between_mins <- matrixStats::rowMins(between_dat, na.rm = T)
    out_dt <- rbind(
      data.table(min_dist = within_mins, comparison = 'within'),
      data.table(min_dist = between_mins, comparison = 'between')
    )[,BC := b]
    return(out_dt)
  })
  spike_dist_dat_out <- rbindlist(spike_dist_dat_out)
  
  p_min_dist <- ggplot2::ggplot(spike_dist_dat_out[min_dist>0], ggplot2::aes(min_dist, fill = comparison)) +
    ggplot2::geom_bar(position = position_dodge()) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Hamming distance to nearest spike-in UMI") +
    ggplot2::ylab("Occurance") +
    ggplot2::theme(legend.position = "top") + 
    ggplot2::geom_vline(xintercept = 2.5, linetype = 'dashed') + 
    ggplot2::xlim(0,15)
  return(p_min_dist)
}