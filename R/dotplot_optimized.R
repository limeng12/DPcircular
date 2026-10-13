#' Dot Plot of Sequence Self-Alignment
#'
#' Computes sliding-window match rates between two sequences and visualizes
#' them as a dot plot, useful for inspecting circular/repeat self-alignments.
#'
#' @param seq1,seq2 Character strings (or DNAString objects coerced to character).
#' @param read_name Axis label for the plot.
#' @param window Sliding window size. Default \code{10}.
#' @param step Step size between windows. Default \code{1}.
#' @param min_match Minimum match rate to keep a point. Default \code{0.9}.
#' @param max_points Maximum number of points to plot (random subsample if exceeded).
#'   Default \code{100000}.
#' @param use_rcpp Logical; use the Rcpp \code{compute_matches} function if available.
#'   Default \code{TRUE}.
#'
#' @return Invisibly returns a \code{data.frame} of plotted points with columns
#'   \code{x}, \code{y}, and \code{match}. Returns \code{NULL} invisibly if no
#'   matches are found.
#'
#' @seealso \code{\link{compute_matches}}, \code{\link{get_circular_count}}
#'
#' @examples
#' \dontrun{
#' dotplot_optimized(dna_seq_char[1], dna_seq_char[1], "read_1")
#' }
#'
#' @export
dotplot_optimized <- function(seq1, seq2, read_name, window = 10, step = 1, min_match = 0.9, 
                              max_points = 100000, use_rcpp = TRUE) {
  
  s1 <- strsplit(as.character(seq1), "")[[1]]
  s2 <- strsplit(as.character(seq2), "")[[1]]
  
  len1 <- length(s1)
  len2 <- length(s2)
  
  # 使用Rcpp如果可用
  if(use_rcpp && require(Rcpp)) {
    # 原来这里 sourceCpp(...) 定义 compute_matches，现改为直接调用包内已编译函数
    match_mat <- compute_matches(s1, s2, window, step, min_match)
    
  } else {
    # 纯R优化版本
    i_pos <- seq(1, len1 - window + 1, step)
    j_pos <- seq(1, len2 - window + 1, step)
    
    match_mat <- matrix(0, nrow = length(i_pos), ncol = length(j_pos))
    
    for(i in seq_along(i_pos)) {
      seq1_window <- s1[i_pos[i]:(i_pos[i] + window - 1)]
      for(j in seq_along(j_pos)) {
        match_mat[i, j] <- sum(seq1_window == s2[j_pos[j]:(j_pos[j] + window - 1)]) / window
      }
    }
    match_mat[match_mat < min_match] <- 0
  }
  
  idx <- which(match_mat > 0, arr.ind = TRUE)
  if(nrow(idx) == 0) return(invisible(NULL))
  
  i_pos <- seq(1, len1 - window + 1, step)
  j_pos <- seq(1, len2 - window + 1, step)
  
  pts <- data.frame(
    x = i_pos[idx[,1]] + window/2,
    y = j_pos[idx[,2]] + window/2,
    match = match_mat[idx]
  )
  
  if(nrow(pts) > max_points) {
    set.seed(42)
    pts <- pts[sample(nrow(pts), max_points), ]
  }
  
  scattermoreplot(pts$x, pts$y, 
                  xlim = c(0, len1), ylim = c(0, len2),
                  cex = 1, col = "black", asp = 1,
                  xlab = read_name, ylab = read_name)
  
  return(invisible(pts))
}



