#' Detect Circular Sequences in a DNA/RNA Sequence
#'
#' Iteratively detects circular (repeat-unit) structures using a modified
#' Smith-Waterman self-alignment. Homopolymer runs >= 6 bp are compressed
#' first; confirmed circles are split at the junction and pushed back for
#' recursive detection. Linear units are aligned with Muscle to build a
#' consensus when more than two are found.
#'
#' @param seq A single DNA/RNA sequence string.
#'
#' @return A named list:
#' \describe{
#'   \item{circle_num}{Number of detected circular units.}
#'   \item{is_circuler}{Logical; TRUE if any circle was found.}
#'   \item{units}{Character vector of linear units.}
#'   \item{consensus}{Consensus of \code{units}, or NULL if <= 2 units.}
#'   \item{score}{Alignment score of the last fragment (only when
#'     \code{consensus} is NULL).}
#' }
#'
#' @seealso \code{\link{sw_circular_cpp}}, \code{\link{filter_false_circular}}
#'
#' @examples
#' \dontrun{
#' res <- get_circular_count(paste0(rep("ATCG", 25), collapse = ""))
#' res$is_circuler
#' }
#'
#' @export
get_circular_count <- function(seq){
  
  # ---- 唯一改动部分：压缩 + 生成freq ----
  chars <- strsplit(seq, "")[[1]]
  r <- rle(chars)
  
  seq_clean_chars <- c()
  freq <- c()
  
  for(i in 1:length(r$lengths)) {
    if(r$lengths[i] >= 6) {
      seq_clean_chars <- c(seq_clean_chars, r$values[i])
      freq <- c(freq, r$lengths[i])
    } else {
      seq_clean_chars <- c(seq_clean_chars, rep(r$values[i], r$lengths[i]))
      freq <- c(freq, rep(1, r$lengths[i]))
    }
  }
  
  seq_clean <- paste0(seq_clean_chars, collapse="")
  # ---- 改动结束 ----
  
  seq_stack <- list(list(clean=seq_clean, freq=freq))
  is_circuler <- FALSE
  circle_num <- 1
  score <- 0
  circular_seqs <- c()
  units <- c()
  
  while(length(seq_stack) > 0) {
    curr <- seq_stack[[1]]
    seq_stack <- seq_stack[-1]
    seq1 <- curr$clean
    freq1 <- curr$freq
    
    if(nchar(seq1) < 15) next
    
    a <- sw_circular_cpp(seq1)
    
    score <- a$score;
    
    
    if(a$is_circular){
      anti_circ <- filter_false_circular(seq1, a$start_j, a$start_i)
      a$is_circular <- a$is_circular && (a$max_score > anti_circ$max_score * 2) # 1.9
    }
    
    if(a$is_circular){
      is_circuler <- TRUE
      circle_num <- circle_num + 1
      split_pos <- a$start_j
      
      # 同步切割
      idx1 <- 1:min(split_pos+5, length(freq1))
      idx2 <- max(1, split_pos-5):min(a$start_i, length(freq1))
      
      aa <- list(clean=substr(seq1, 1, max(idx1)), freq=freq1[idx1])
      bb <- list(clean=substr(seq1, min(idx2), max(idx2)), freq=freq1[idx2])
      
      seq_stack <- c(seq_stack, list(aa, bb))
      circular_seqs <- c(circular_seqs, paste0(seq1, "$", aa$clean, "$", bb$clean, "$", circle_num))
      
    } else {
      # 还原
      units <- c(units, paste0(sapply(1:length(freq1), function(i) 
        paste0(rep(substr(seq1,i,i), freq1[i]), collapse="")), collapse=""))
    }
  }
  
  #  circular_seqs=circular_seqs,
  
  #3个untis以上才有必要求取consensus
  if(length(units) > 2){
    
    alignment <- msa( DNAStringSet(units), method = "Muscle")
    consensus <- consensusString(alignment, ambiguityMap=IUPAC_CODE_MAP)
    consensus <-  str_replace_all(consensus, "-", "")
    
    return(list(circle_num=circle_num, is_circuler=is_circuler,
                units=units, consensus=consensus) )
  
  }else{
    return(list(circle_num=circle_num, score = score, 
                is_circuler=is_circuler, 
                units=units, consensus=NULL) )
    
  }
  
  
}




filter_false_circular <- function(seq_clean, start_j, end_pos) {
  
  d <- str_sub(seq_clean, 1, start_j)
  e <- str_sub(seq_clean,1, end_pos)
  d_rev <- stri_reverse(d)
  dd<-sw_cpp(d_rev, e);
  
  return(dd)
}

