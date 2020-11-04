#' Load chromosome length and bin counts into global environment
#'
#' @param chrlen.file Text file for chromosome length
#' @param bin.width The bin width.
#' @export
load.chrlen <- function(chrlen_vec, bin.width) {
        if (exists("chr.len", envir = .GlobalEnv)) {
                if (identical(chr.len, chrlen_df$X2)) {
                        if (ceiling(chr.len[1] / bin.counts[1]) == bin.width)
                                return(NULL)
                }
        }
        # when given a file
        # scan(chrlen.file, quiet = T)
        message(paste(
                "loading chromosome bin lengths for bin width:",
                bin.width
        ))
        assign("chr.len", chrlen_vec, envir = .GlobalEnv)
        assign("bin.counts", ceiling(chr.len / bin.width), envir = .GlobalEnv)
        assign("bin.from", c(0, cumsum(bin.counts)), envir = .GlobalEnv)
}

#' Create a \code{\link{GRanges}} object for bins.
#'
#' @param chrlen_df data.frame for chromosome sizes
#' @param bin_width Bin width
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
make_gr_bins <- function (chrlen_df, bin_width) {
        load.chrlen(chrlen_df$X2, bin_width)
        range.start <- unlist(lapply(bin.counts[1:nrow(chrlen_df)],
                                     function(count) {
                                             seq(from = 1,
                                                 by = bin_width,
                                                 length = count)
                                     }))
        range.end <- range.start + bin_width - 1
        for (chr in 1:nrow(chrlen_df)) {
                range.end[bin.from[chr + 1]] <- chrlen_df$X2[chr]
        }
        gr <- GRanges(
                seqnames = Rle(chrlen_df$X1,
                               bin.counts[1:nrow(chrlen_df)]),
                ranges = IRanges(start = range.start,
                                 end = range.end)
        )
        return(gr)
}

#' Convert base pair indices to bin indices
#'
#' @param bin Base pair indices
#' @param bin.width Bin width
#' @export
bp2bin <- function(bp, bin.width) {
        ceiling(bp / bin.width)
}

#' Convert bin indices to base pair indices
#'
#' @param bin Bin indices
#' @param bin.width Bin width
#' @export
bin2bp <- function(bin, bin.width = 50) {
        (bin - 1) * bin.width + 1
}

#' Imports narrowPeak files
#'
#' @param f File path
#' @export
#' @importFrom rtracklayer import.bed
import.narrowPeak <- function(f) {
        extraCols.narrowPeak <-
                c(
                        singnalValue = "numeric",
                        pValue = "numeric",
                        qValue = "numeric",
                        peak = "integer"
                )
        import.bed(f, extraCols = extraCols.narrowPeak)
}

#' Imports broadPeak files
#'
#' @param f File path
#' @export
#' @importFrom rtracklayer import.bed
import.broadPeak <- function(f) {
        extraCols_broadPeak <- c(
                signalValue = "numeric",
                pValue = "numeric",
                qValue = "numeric"
                )
        import.bed(f, extraCols = extraCols_broadPeak)
}

#' Convert seqnames to numeric chromosome indices. Old. For the new one, see seq2num.
#'
#' @export
chr2num <- function(c) {
        # earlier version need to merge with seq2num
        # does not handle Y M yet
        what <- gsub("chr", "", as.character(c))
        if (what == "X") {
                return(23)
        } else {
                return(as.numeric(what))
        }
}

#' Convert seqnames to numeric chromosome indices. chrX mapped to 23; chrM and chrY mapped to NAs.
#'
#' @export
seq2num <- function(seq, chrlen_df) {
        if (is.character(seq))
                seq <- factor(seq, levels = chrlen_df$X1)
        as.numeric(seq)
}

#' Convert vector of genome wide bin indices to a list of indices per chromosome.
#'
#' @param bins Vector of bin indices
#' @param bin.width Bin width
#' @export
gw2chr <- function(bins, bin.width = 200) {
        # input: genome wide incices of bins
        # output: list with an item for each chromosome
        if (!exists("bin.from") | !exists("bin.counts")) {
                stop("please load chromosome lengths first")
        }
        out <- list()
        for (chr in 1:23) {
                out[[chr]] <-
                        bins[(bins <= bin.from[chr + 1]) &
                                     (bins > bin.from[chr])] -
                        bin.from[chr]
        }
        return(out)
}

#' Convert a list of indices per chromosome (must have 23 items in list) to a vector of genome wide bin indices
#'
#' @param bins Vector of bin indices
#' @param bin.width Bin width
#' @export
list2gw <- function(bins.ls, bin.width = 200) {
        # input: list of 23 chromosomes
        # output: output a long vector of bin indices
        if (!exists("bin.from") | !exists("bin.counts")) {
                stop("please load chr.len first")
        }
        if (length(bins.ls) != 23) {
                stop("list length is not 23")
        }
        unlist(lapply(1:23, function(chr)
                bins.ls[[chr]] + bin.from[chr]))
}

#' Convert chromosome bin indices to genome wide bin indices.
#'
#' @param chr vector of chr (Rle)
#' @param bins vector of bins
#' @export
chr2gw <- function(chr, bins, chrlen_df, bin.width) {
        # check if chr is numeric
        load.chrlen(chrlen_df$X2, bin.width)
        if (!is.numeric(chr))
                chr <- seq2num(chr, chrlen_df)
        bins <- bins[!is.na(chr)]
        chr <- chr[!is.na(chr)]
        bin.from[chr] + bins
}

#' Average matrix rows based on a factor
#'
#' @export
aveMatFac <- function(mat, fac) {
        # need to be able to handle character or numeric
        if (class(fac) != "factor")
                fac <- factor(fac)
        rown <- length(levels(fac))
        coln <- dim(mat)[2]
        out <- matrix(, rown, coln)
        ind <- as.numeric(fac)
        for (i in 1:rown) {
                out[i,] <- colMeans(mat[ind == i, , drop = F], na.rm = T)
        }
        rownames(out) <- levels(fac)
        return(out)
}

#' Sum matrix rows based on a factor
#'
#' @export
sumMatFac <- function(mat, fac) {
        # need to be able to handle character or numeric
        if (class(fac) != "factor")
                fac <- factor(fac)
        rown <- length(levels(fac))
        coln <- dim(mat)[2]
        out <- matrix(, rown, coln)
        ind <- as.numeric(fac)
        for (i in 1:rown) {
                out[i,] <- colSums(mat[ind == i, , drop = F], na.rm = T)
        }
        rownames(out) <- levels(fac)
        return(out)
}

#' Permutate columns of matrix independently
#'
#' @export
perm.mat <- function(mat) {
        out <- matrix(, nrow(mat), ncol(mat))
        for (i in 1:ncol(mat)) {
                out[, i] <- sample(mat[, i])
        }
        return(out)
}

#' Permutate columns of multiple matrix independently
#'
#' @export
perm.mat.multi <- function(...) {
        mat.list <- list(...)
        # to-do: make sure elements are all matrices and have then same dimension
        n.col <- ncol(mat.list[[1]])
        n.row <- nrow(mat.list[[2]])
        out <- list()
        for (j in 1:length(mat.list)) {
                out[[j]] <- matrix(, nrow(mat.list[[j]]), ncol(mat.list[[j]]))
        }
        for (i in 1:n.col) {
                pidx <- sample(1:n.row)
                for (j in 1:length(mat.list)) {
                        out[[j]][, i] <- mat.list[[j]][pidx, i]
                }
                
        }
        return(out)
}

#' Count reads in bins from GeonomicAlignments object.
#'
#' @param align A \code{\link{GenomicAlignments}} object.
#' @param paried A logical value indicating if alignment object is paired.
#' @export
#' @importFrom GenomicRanges seqnames start end
countReads <-
        function(align,
                 chrlen.file,
                 chr.count,
                 bin.width,
                 paired = F,
                 counts = NULL) {
                load.chrlen(chrlen.file, bin.width)
                if (is.null(counts))
                        counts <- numeric(bin.from[chr.count + 1])
                
                if (paired) {
                        rstart <- pmin(start(align@first), start(align@last))
                        rend <- pmax(end(align@first), end(align@last))
                } else {
                        rstart <- start(align)
                        rend <- end(align)
                }
                ctr <- ceiling((rstart + rend) / 2)
                bin.gw <- chr2gw(seqnames(align), bp2bin(ctr, bin.width))
                
                count_bins(counts, bin.gw)
        }

#' Convert bam files to bin counts
#'
#' @param bam.file The bam file to be converted.
#' @export
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @importFrom Rsamtools BamFile
bam2bin <-
        function(bam.file,
                 chrlen.file,
                 chr.count,
                 bin.width,
                 paired = F,
                 ...) {
                if (is.character(bam.file)) {
                        bam.file <- BamFile(bam.file, ...)
                }
                open(bam.file)
                out <- NULL
                repeat {
                        if (paired) {
                                alignment <- readGAlignmentPairs(bam.file)
                        } else {
                                alignment <- readGAlignments(bam.file)
                        }
                        if (length(alignment) == 0L) {
                                break
                        }
                        out <-
                                countReads(
                                        alignment,
                                        chrlen.file = chrlen.file,
                                        chr.count = chr.count,
                                        bin.width = bin.width,
                                        paired = paired,
                                        counts = out
                                )
                }
                close(bam.file)
                out
        }

#' Convert reads saved as list of data.frames to counts
#'
#' @export
df2bin <-
        function(df,
                 chrlen_df,
                 bin_width,
                 paired = F,
                 counts = NULL) {
                stopifnot(all(chrlen_df$X1 %in% names(df)))
                df <- df[chrlen_df$X1]
                
                load.chrlen(chrlen_df$X2, bin_width)
                reads_gw <- lapply(1:nrow(chrlen_df), function(i) {
                        if (paired) {
                                stopifnot(
                                        names(df[[i]]) == c(
                                                "firststart",
                                                "firstend",
                                                "laststart",
                                                "lastend"
                                        )
                                )
                                rstart <- pmin(df[[i]]$firststart, df[[i]]$laststart)
                                rend <- pmax(df[[i]]$firstend, df[[i]]$lastend)
                        } else {
                                stopifnot(names(df[[i]]) == c("start", "end"))
                                rstart <- df[[i]]$start
                                rend <- df[[i]]$end
                        }
                        
                        chr2gw(chrlen_df$X1[i],
                               bp2bin(ceiling((
                                       rstart + rend
                               ) / 2), bin_width),
                               chrlen_df,
                               bin_width)
                })
                reads_gw_vec <- unlist(reads_gw)
                
                if (is.null(counts))
                        counts <- numeric(bin.from[nrow(chrlen_df) + 1])
                count_bins(counts, reads_gw_vec)
        }

#' Convert tagAlign files to bin counts
#'
#' @param ta.file bam The tagAlign file to be converted.
#' @export
#' @importFrom rtracklayer import.bed
ta2bin <- function(ta.file, ...) {
        alignment <- import.bed(ta.file)
        message(paste("number of tags:", length(alignment)))
        countReads(alignment, ...)
}
#' Convert string to GR object
#'
#' @export
str_to_gr <- function(regions, sep = c(":", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- tidyr::separate(
    data = ranges.df,
    col = 'ranges',
    sep = paste0(sep[1], "|", sep[2]),
    into = c('chr', 'start', 'end')
  )
  granges <- GenomicRanges::makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

#' Count overlaps of regions and read centers
#' @export
df2region <- function(df, chrlen_df, regions_gr, paired = F) {
    stopifnot(all(chrlen_df$X1 %in% names(df)))
    reads_gr <- do.call(c, lapply(1:nrow(chrlen_df), function(i) {
            if (paired) {
                    stopifnot(
                            names(df[[i]]) == c(
                                    "firststart",
                                    "firstend",
                                    "laststart",
                                    "lastend"
                            )
                    )
                    rstart <- pmin(df[[i]]$firststart, df[[i]]$laststart)
                    rend <- pmax(df[[i]]$firstend, df[[i]]$lastend)
            } else {
                    stopifnot(names(df[[i]]) == c("start", "end"))
                    rstart <- df[[i]]$start
                    rend <- df[[i]]$end
            }
            reads_ctr = ceiling((
                                   rstart + rend
                                                  ) / 2)
            gr_out = GenomicRanges::GRanges(seqnames = chrlen_df$X1[i],
        IRanges::IRanges(start = reads_ctr,
            end = reads_ctr))
            gr_out
            }))
    counts_out = GenomicRanges::countOverlaps(regions_gr, reads_gr)
    counts_out
}

#' Calculate empirical p-values
#'
#' @param null.vec Vector of empirical null distribution
#' @param obs.vec Vector of observed test statistics
#' @param alternative Type of alternative hypothesis
#' @export
empirical.pvalue <-
        function(null.vec,
                 obs.vec,
                 alternative = c("less", "greater", "two.sided")) {
                # given a vector of null distribution, and a vector of observed statistics
                # calculate a vector of p-values
                # allows two sided p-values
                # how to handle NAs in null?
                
                n <- length(obs.vec)
                pvalue.vec <- numeric(n)
                alternative = match.arg(alternative)
                null.med <- median(null.vec, na.rm = T)
                for (i in 1:n) {
                        if (is.na(obs.vec[i])) {
                                pvalue.vec[i] <- NA
                        } else {
                                pvalue.vec[i] <-
                                        switch(
                                                alternative,
                                                "less" = {
                                                        (sum(null.vec <= obs.vec[i], na.rm = T) + 1) /
                                                                (length(null.vec) + 1)
                                                },
                                                "greater" = {
                                                        (sum(null.vec >= obs.vec[i], na.rm = T) + 1) /
                                                                (length(null.vec) + 1)
                                                },
                                                "two.sided" = {
                                                        ifelse(
                                                                obs.vec[i] > null.med,
                                                                2 * (
                                                                        sum(null.vec >= obs.vec[i], na.rm = T) + 1
                                                                ) /
                                                                        (length(
                                                                                null.vec
                                                                        ) + 1),
                                                                2 * (
                                                                        sum(null.vec <= obs.vec[i], na.rm = T) + 1
                                                                ) /
                                                                        (length(
                                                                                null.vec
                                                                        ) + 1)
                                                        )
                                                }
                                        )
                                pvalue.vec[i] <- pmin(pvalue.vec[1], 1)
                        }
                }
                return(pvalue.vec)
        }