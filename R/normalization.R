#' normalized_expression_data
#'
#' A function to compute the mean of a vector
#' @param expression_data Expression matrix: cells * genes
#' @export
#' @examples
#' mean(1:3)
#' \dontrun{ mean(1:1e99) }
#'

normalized_expression_data <- function(expression_data) {
  normalized_counts <- apply(expression_data, 2, function(x) {
    if (all(x == 0)) {
      return(x)
    } else if (all(x == x[1])) {
      return(rep(1, length(x)))
    } else if (all(x >= 0)) {
      min_val <- min(x[x > 0])
      return(ifelse(x == 0, 0, (x - min_val) / (max(x) - min_val)))
    } else {
      return((x - min(x)) / (max(x) - min(x)))
    }
  })
  normalized_counts <- normalized_counts
  return(normalized_counts)
}
