#' @title A function converting a portfolio-stock dataframe to a stock-count matrix
#' @description A function converting a portfolio-stock dataframe to a stock-count matrix
#' @param d data frame with each row containing the first 10 core stocks of a portfolio
#' @return a list containing the converted matrix and whether there are false data
#' @export
convert <- function(d){
  t <- d
  rownames(t) <- NULL
  t <- t[!duplicated(t), ]
  C <- dim(t)
  vec <- numeric(C[1] * C[2])
  for(i in 1:C[1]){
    for(j in 1:C[2]){
      vec[(i - 1) * C[2] + j] <- t[i, j]
    }
  }
  stock <- unique(vec)
  m <- matrix(0, nrow = C[1], ncol = length(stock))
  for(i in 1:nrow(m)){
    for(j in 1:C[2]){
      m[i, ] <- m[i, ] + (stock == t[i, j])
    }
  }
  colnames(m) <- stock
  count1 <- sum(m == 1)
  count0 <- sum(m == 0)
  m1 <- m
  stockorder <- names(sort(colSums(m1), decreasing = T))
  m2 <- matrix(0, nrow = C[1], ncol = length(stockorder))
  for(i in 1:nrow(m2)){
    for(j in 1:C[2]){
      m2[i, ] <- m2[i, ] + (stockorder == t[i, j])
    }
  }
  colnames(m2) <- stockorder
  flag <- ((count1+count0) == nrow(m)*ncol(m))
  result <- list(res = m2, data = t, mat.m = m, correct = flag)
  return(result)
}

#' @title correct errors from the result of the function convert
#' @description if there are false data in the result it, find it and correct it
#' @param list result of the function convert
#' @return return a correct matrix 
#' @export
correct.error <- function(list){
  m <- list$mat.m
  err <- which(m == 2) %% nrow(m)
  t_fix <- list$data[-c(err), ]
  D <- dim(t_fix)
  vec_star <- numeric(D[1] * D[2])
  for(i in 1:D[1]){
    for(j in 1:D[2]){
      vec_star[(i - 1)*D[2] + j] <- t_fix[i, j]
    }
  }
  stock_star <- unique(vec_star)
  m3 <- matrix(0, nrow = D[1], ncol = length(stock_star))
  for(i in 1:nrow(m3)){
    for(j in 1:D[2]){
      m3[i, ] <- m3[i, ] + (stock_star == t_fix[i, j])
    }
  }
  colnames(m3) <- stock_star
  count1 <- sum(m3 == 1)
  count0 <- sum(m3 == 0)
  stock_star_order <- names(sort(colSums(m3), decreasing = T))
  m4 <- matrix(0, nrow = D[1], ncol = length(stock_star_order))
  for(i in 1:nrow(m4)){
    for(j in 1:D[2]){
      m4[i, ] <- m4[i, ] + (stock_star_order == t_fix[i, j])
    }
  }
  colnames(m4) <- stock_star_order
  flag <- ((count1 + count0) == nrow(m3) * ncol(m3))
  result <- list(res = m4, correct = flag)
  return(result)
}

#' @title A dataset used to describe how to use the function convert
#' @name lmlv19
#' @description This dataset contains no false data, so it can be directly converted by function convert
#' @examples
#' \dontrun{
#' data(lmlv19)
#' attach(lmlv19)
#' l <- convert(lmlv19)
#' str(l)
#' }
NULL

#' @title A dataset used to describe how to use the function correct.error
#' @name smsv19
#' @description This dataset contains errors, so it can be corrected by function correct.error
#' @examples
#' \dontrun{
#' data(smsv19)
#' attach(smsv19)
#' l <- convert(smsv19)
#' res <- correct.error(l)
#' str(res)
#' }
#' @import bootstrap
#' @import stats4
#' @import boot
#' @importFrom Rcpp evalCpp
#' @useDynLib SA23204171
NULL

#' @title The dataset ironslag
#' @name ironslag
#' @description
#' See the data for details
NULL