## -----------------------------------------------------------------------------
# stock name
stock <- LETTERS
l <- list()
for(i in 1:10){
  # generate 10 portfolios
  l[[i]] <- sample(stock, 10)
}
# portfolio is a data frame with each row containing the 10 core stocks of a 
# portfolio. portfolio has no error data, because every stock appears at most 
# once in each row.
portfolio <- as.data.frame(rbind(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]],
                                 l[[7]], l[[8]], l[[9]], l[[10]]))
# convert the data frame
library(SA23204171)
result <- convert(portfolio)

## -----------------------------------------------------------------------------
# stock names
stock <- LETTERS
l <- list()
for(i in 1:10){
  # generate 10 portfolios
  l[[i]] <- sample(stock, 10)
}
# this portfolio contains false data, with a few stocks appearing more than once
# in the same portfolio
l[[1]][10] <- l[[1]][5]
l[[6]][8] <- l[[6]][2]
portfolio <- as.data.frame(rbind(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]],
                                 l[[7]], l[[8]], l[[9]], l[[10]]))
# first use the function convert() to convert the data frame
library(SA23204171)
result <- convert(portfolio)
result$correct
# FALSE means there are errors in the original data, so we use correct.error() 
# to filter out these false portfolios and give the correct converted matrix.
# The error is apparent because there are "2" in the converted matrix.
new_result <- correct.error(result)
new_result

