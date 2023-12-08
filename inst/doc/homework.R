## -----------------------------------------------------------------------------
x <- 1
y <- 2
print(x+y)

## ----eval=FALSE---------------------------------------------------------------
#  x <- 1
#  y <- 2
#  print(x+y)

## ----echo=FALSE---------------------------------------------------------------
x <- 1
y <- 2
print(x+y)

## -----------------------------------------------------------------------------
my.sample <- function(x, n, prob = NULL){
  ### 给定一有限总体，有放回地从中抽取样本
  #   x: 输入一个有限总体，为向量形式
  #   n: 从x中产生的随机数个数，为一正整数
  #   prob: 与x等长的向量，表示x中对应位次元素的入样概率，未指
  #定则默认等概率抽取
  
  if(is.null(prob)){
    # 未指定prob，默认等概率抽取
    p <- 1:length(x)/sum(1:length(x))
    cp <- cumsum(p)
    U <- runif(n)
    r <- x[findInterval(U, cp)+1]
  } else{
    # 用户指定了各元素入样概率
      if((length(prob) != length(x)) | (sum(prob) != 1)){
      # 用户指定的prob格式错误，则返回报错信息
      return("Wrong entry of the argument 'prob'!")
    } else{
      p <- cumsum(prob)
      U <- runif(n)
      r <- x[findInterval(U, p)+1]
    }
  }
  
  return(r) # 返回抽样结果，为一向量
}

# tests
x <- c(1, 2, 4, 10)
n <- 10
p1 <- 1:4/10
p2 <- c(0.3, 0.3, 0.4) # 错误指定prob
p3 <- c(0.1, 0.1, 0.3, 0.2, 0.3) # 错误指定prob
my.sample(x, n, p1)
my.sample(x, n, p2)
my.sample(x, n, p3)

# comparison with 'sample()'
x1 <- 1:10
pr <- 1:10/sum(1:10)
N <- 10000
res1 <- sample(x1, N, prob = pr, replace = T)
res2 <- my.sample(x1, N, prob = pr)
counts <- table(c(rep(0, N), rep(1, N)), c(res1, res2))
barplot(counts, 
        main="my.sample() vs sample()", 
        ylab='Count', 
        col=c("blue","red"), 
        beside=TRUE)
legend('topleft', 
       fill = c('blue', 'red'), 
       col = c('blue', 'red'), 
       legend = c('sample', 'my.sample'))

## -----------------------------------------------------------------------------
u <- runif(1000)
x <- ifelse(u <= 0.5, log(2*u), -log(2-2*u))
# 对生成的随机数做直方图
hist(x, 
     freq = F, 
     xlim = c(-10, 10), 
     ylim = c(0, 0.5))
# 叠加标准Laplace分布的密度曲线
curve(0.5*exp(-abs(x)), 
      col = 'red', 
      lwd = 1.5, 
      from = -10, to = 10, 
      add = T)

## -----------------------------------------------------------------------------
# generation function
my.rbeta <- function(a, b, size){
  
  n <- size
  result <- numeric(n) # 生成的随机序列储存在向量result中
  flag <- 1 # 计数变量，一旦大于n则退出whlie循环
  
  # rho()函数
  rho <- function(x){
    m1 <- (a+b-2)^(a+b-2) * x^(a-1) * (1-x)^(b-1)
    m2 <- (a-1)^(a-1) * (b-1)^(b-1)
    return(m1/m2)
  }
  
  while(flag <= n){
    u <- runif(1) # 生成一个U(0, 1)上的随机数
    y <- runif(1) # 生成来自包络分布的一个随机数
    if(u <= rho(y)){
      result[flag] <- y
      flag <- flag + 1
    }
  }
  
  return(result)
}

# tests
res <- my.rbeta(3, 2, 1000)
# 对生成的随机数做直方图
hist(res, 
     probability = T , 
     xlim = c(0, 1))
# 叠加Beta(3, 2)分布的密度曲线
curve(12* x^2 * (1-x), 
      col = 'red', 
      lwd = 1.5, 
      from = 0, to = 1, 
      add = T)

## -----------------------------------------------------------------------------
# generation function
my.repanechnikov <- function(size){
  
  n <- size
  result <- numeric(n) # 生成的随机序列储存在向量result中
  
  for(i in 1:n){
    r <- runif(3, -1, 1)
    if((abs(r[3]) >= abs(r[2])) & (abs(r[3]) >= abs(r[1]))){
      result[i] <- r[2]
    } else{
      result[i] <- r[3]
    }
  }
  
  return(result)
}

# tests
n <- 10000
res <- my.repanechnikov(n)
# 对生成的随机数做直方图
hist(res, 
     probability = T , 
     xlim = c(-1, 1))
# 叠加Epanechnikov核密度曲线
curve(3*(1-x^2)/4, 
      col = 'red', 
      lwd = 1.5, 
      from = -1, to = 1, 
      add = T)

## -----------------------------------------------------------------------------
# rho_min = 1, 另外两个rho为0.5、0.8
rho <- c(1, 0.5, 0.8)
# 固定两平行线之间的距离d, 针长l = rho * d的改变通过rho的改变来实现
d <- 1
# 重复模拟次数
K <- 100 
# 每次模拟的投针次数
n <- 1e6
# 模拟结果存储在该列表中
result <- list("rho_is_1" = numeric(K), 
               "rho_is_0.5" = numeric(K), 
               "rho_is_0.8" = numeric(K))

for(i in 1:3){
  for(k in 1:K){
    l <- d*rho[i]
    X <- runif(n, 0, d/2)
    Y <- runif(n, 0, pi/2)
    pi_hat <- 2*rho[i]/mean((l/2*sin(Y))>X)
    result[[i]][k] <- pi_hat
  }
}

# rho = 1的方差
var(result$rho_is_1)
# rho = 0.5的方差
var(result$rho_is_0.5)
# rho = 0.8的方差
var(result$rho_is_0.8)

## -----------------------------------------------------------------------------
# 产生N个均匀随机变量用以计算MC估计量
N <- 2000
# 重复模拟K次
K <- 500

# simple MC method
theta1_hat <- numeric(K)
for(i in 1:K){
  u <- runif(N)
  theta1_hat[i] <- mean(exp(u))
}
va1 <- var(theta1_hat); va1

# antithetic variate approach
theta2_hat <- numeric(K)
for(i in 1:K){
  u <- runif(N/2)
  theta2_hat[i] <- sum(exp(u) + exp(1-u))/N
}
va2 <- var(theta2_hat); va2

# 方差减少的百分比
cat(round((va1-va2)*100/va1, 2), '%', sep = '')

## -----------------------------------------------------------------------------
N=100 #重复模拟100次
m=1000 #每次产生的随机数个数
res1 <- numeric(N) #用f_1估计的结果向量
res2 <- numeric(N) #用f_2估计的结果向量

# 3个函数
f1 <- function(x) exp(1-x) * (x>1)
f2 <- function(x) x * exp((1-x^2)/2) * (x>1)
g <- function(x) x^2 * exp(-x^2/2) * (x>1) / sqrt(2*pi)

# 用f_1进行估计
for (i in 1:N){
  # 用逆变换法产生m个来自f1分布的随机数
  f1_rev <- function(x) 1 - log(1-x)
  u <- runif(m)
  x <- f1_rev(u)
  res1[i] <- sum(g(x)/f1(x))/m
}

# 用f_2进行估计
for (i in 1:N){
  # 用逆变换法产生m个来自f2分布的随机数
  f2_rev <- function(x) sqrt(1 - 2*log(1-x))
  u <- runif(m)
  x <- f2_rev(u)
  res2[i] <- sum(g(x)/f2(x))/m
}

# 求方差, 可以看出用f2求出的估计方差更小
var(res1)
var(res2)

# 作图分别查看f1与g、f2与g的正比关系:
f1_g <- function(x) f1(x)/g(x)
f2_g <- function(x) f2(x)/g(x)
curve(f1_g(x), 
      from = 1, to = 5, 
      xlab = 'x', ylab = 'proportion')
curve(f2_g(x), 
      from = 1, to = 5, 
      col = 'red', 
      add = T)
legend('topleft',
       lty = c(1, 1), 
       col = c('black', 'red'), 
       legend = c('f1 vs g', 'f2 vs g'))
# 显然f2与g的正比关系更强

## -----------------------------------------------------------------------------
head(res1)
head(res2)

## -----------------------------------------------------------------------------
# 分成k = 5个区间
k <- 5
# 每个区间生成5000个随机数
m <- 1000
# 目标参数的被积函数
g <- function(x) (exp(-x))/(1+x^2)
# 第j个区间的条件密度函数f_j
fj <- function(x, j){
  5 * exp(j-1-5*x)*(x > ((j-1)/5))*(x < (j/5)) / (1-exp(-1))
}
# f_j对应的分布函数
Fj <- function(x, j) (1 - exp(j-1-5*x)) / (1 - exp(-1))
# f_j对应的分布函数的反函数
Fj_rev <- function(x, j) (j - log(exp(1) - (exp(1)-1)*x))/5

# 初始化
theta_hat <- 0

# stratified importance sampling estimation
for(j in 1:k){
  # 逆变换法生成来自密度f_j的随机数
  u <- runif(m)
  x <- Fj_rev(u, j)
  # 逐个区间进行分段积分的估计, 然后累加得到全区间积分的最终估计
  theta_hat <- theta_hat + sum(g(x)/fj(x, j))/m
}

# 最终估计值, 可以看出与Example 5.10的结果是非常接近的
theta_hat

## -----------------------------------------------------------------------------
# 重复模拟m次
m <- 1e4
# 每次模拟的样本容量为n
n <- 20
# 自由度为n-1的t分布的上0.025分位数(均值采用双边检验)
c1 <- qt(0.975, n-1)
# 自由度为n-1的卡方分布的下0.05分位数(方差采用单边检验)
c2 <- qchisq(0.05, n-1)
# 指示向量，用来指示每次置信区间是否包含目标参数
res_mean <- numeric(m)
res_var <- numeric(m)

for(i in 1:m){
  # 生成简单卡方样本
  x <- rchisq(n, 2)
  
  # 卡方均值的t区间估计
  mean_lower_hat <- mean(x) - c1 * sqrt(var(x)) / sqrt(n)
  mean_upper_hat <- mean(x) + c1 * sqrt(var(x)) / sqrt(n)
  # 本次t区间估计是否覆盖真实卡方均值
  res_mean[i] <- ifelse((mean_lower_hat <= 2)&(mean_upper_hat >= 2), 1, 0)
  
  # 卡方方差的区间估计
  var_lower_hat <- (n-1)*var(x)/c2
  # 本次区间估计是否覆盖真实卡方方差
  res_var[i] <- ifelse(var_lower_hat >= 4, 1, 0)
}
# 均值的覆盖概率估计
mean_CP_hat <- mean(res_mean)
mean_CP_hat
# 方差的覆盖概率估计
var_CP_hat <- mean(res_var)
var_CP_hat

## -----------------------------------------------------------------------------
M <- 1000 # 1000次重复模拟
a <- 0.1 # 显著性水平alpha
b1 <- 0.1; b2 <- 1 # beta分布的参数
p_bon <- list() # 经过Bonferroni校正的1000组p值结果列表
p_BH <- list() # 经过BH校正的1000组p值结果列表

# 1000次校正过程
for (i in 1:M){
  p <- c(runif(950), rbeta(50, b1, b2))
  p_bon[[i]] <- p.adjust(p, method = 'bonferroni')
  p_BH[[i]] <- p.adjust(p, method = 'BH')
}

# 两种校正后的FWER估计值
FWER_bon_vec <- numeric(M)
FWER_BH_vec <- numeric(M)
for (i in 1:M){
  FWER_bon_vec[i] <- (sum(p_bon[[i]][1:950] < a) > 0)
  FWER_BH_vec[i] <- (sum(p_BH[[i]][1:950] < a) > 0)
}
FWER_bon <- mean(FWER_bon_vec)
FWER_BH <- mean(FWER_BH_vec)

# 两种校正后的FDR估计值
FDR_bon_vec <- numeric(M)
FDR_BH_vec <- numeric(M)
for (i in 1:M){
  FDR_bon_vec[i] <- (sum(p_bon[[i]][1:950] < a)/sum(p_bon[[i]] < a))
  FDR_BH_vec[i] <- (sum(p_BH[[i]][1:950] < a)/sum(p_BH[[i]] < a))
}
FDR_bon <- mean(FDR_bon_vec)
FDR_BH <- mean(FDR_BH_vec)

# 两种校正后的TPR估计值
TPR_bon_vec <- numeric(M)
TPR_BH_vec <- numeric(M)
for (i in 1:M){
  TPR_bon_vec[i] <- (sum(p_bon[[i]][951:1000] < a)/50)
  TPR_BH_vec[i] <- (sum(p_BH[[i]][951:1000] < a)/50)
}
TPR_bon <- mean(TPR_bon_vec)
TPR_BH <- mean(TPR_BH_vec)

# results of Bonferroni adjustment
print(round(c(FWER_bon, FDR_bon, TPR_bon), 3))
# results of B-H adjustment
print(round(c(FWER_BH, FDR_BH, TPR_BH), 3))

## -----------------------------------------------------------------------------
m <- 1000 # 重复模拟1000次
B <- 1000 # 每次模拟bootstrap采样次数
n <- c(5, 10, 20) # 样本量大小

# 各参数理论值
lambda <- 2
bias <- lambda/(n-1)
se <- lambda*n/((n-1)*sqrt(n-2))

# bootstrap procedure

# 3种样本量
x1 <- rexp(n[1], lambda)
x2 <- rexp(n[2], lambda)
x3 <- rexp(n[3], lambda)
# 参数lambda在各样本下的m次bootstrap估计
lambda1_boot <- numeric(m)
lambda2_boot <- numeric(m)
lambda3_boot <- numeric(m)
# 偏差在各样本下的m次bootstrap估计
bias1_boot <- numeric(m)
bias2_boot <- numeric(m)
bias3_boot <- numeric(m)
# 标准差在各样本下的m次bootstrap估计
se1_boot <- numeric(m)
se2_boot <- numeric(m)
se3_boot <- numeric(m)

for(i in 1:m){
  # 第i次重复模拟时各样本的lambda在B次重采样中的bootstrap估计值
  lambda1_temp <- numeric(B)
  lambda2_temp <- numeric(B)
  lambda3_temp <- numeric(B)
  
  # 进行B次重采样
  for(j in 1:B){
    # 对3个样本分别生成重采样样本
    x1_boot <- sample(x1, replace = T)
    x2_boot <- sample(x2, replace = T)
    x3_boot <- sample(x3, replace = T)
    # 3个样本在第j次重采样的lambda估计
    lambda1_temp[j] <- 1/mean(x1_boot)
    lambda2_temp[j] <- 1/mean(x2_boot)
    lambda3_temp[j] <- 1/mean(x3_boot)
  }
  
  # 取B次重采样lambda估计的均值作为lambda的第i次bootstrap估计
  lambda1_boot[i] <- mean(lambda1_temp)
  lambda2_boot[i] <- mean(lambda2_temp)
  lambda3_boot[i] <- mean(lambda3_temp)
  # 取B次lambda估计的均值与原样本lambda估计的差作为bias的第i次bootstrap估计
  bias1_boot[i] <- lambda1_boot[i] - 1/mean(x1)
  bias2_boot[i] <- lambda2_boot[i] - 1/mean(x2)
  bias3_boot[i] <- lambda3_boot[i] - 1/mean(x3)
  # 取B次lambda估计的标准差作为se的第i次bootstrap估计
  se1_boot[i] <- sd(lambda1_temp)
  se2_boot[i] <- sd(lambda2_temp)
  se3_boot[i] <- sd(lambda3_temp)
}

# comparison between mean bootstrap bias and the theoretical one
print(bias)
print(c(mean(bias1_boot), mean(bias2_boot), mean(bias3_boot)))
# comparison between mean bootstrap standard error and the theoretical one
print(se)
print(c(mean(se1_boot), mean(se2_boot), mean(se3_boot)))

## -----------------------------------------------------------------------------
library(bootstrap)
str(law)

B <- 500 # 进行B次重采样(一阶、二阶)
R1 <- numeric(B) # 一阶重采样中相关系数的bootstrap估计
se_R1 <- numeric(B) # 一阶重采样相关系数bootstrap估计的标准差估计
R_hat <- cor(law$LSAT, law$GPA) # 原始样本的相关系数估计
alpha <- 0.05 # 置信区间置信度为95%

# bootstrap
for(i in 1:B){
  # 第i次重采样样本
  id1 <- sample(1:nrow(law), replace = T)
  # 第i次重采样中相关系数的bootstrap估计
  R1[i] <- cor(law$LSAT[id1], law$GPA[id1])
  
  # 下面需要利用bootstrap估计第i次重采样样本相关系数估计的标准差
  
  # 对第i次重采样样本进行二阶重采样
  R2 <- numeric(B) # 二阶重采样相关系数的bootstrap估计
  for(j in 1:B){
    id2 <- sample(id1, replace = T)
    R2[j] <- cor(law$LSAT[id2], law$GPA[id2])
  }
  se_R1[i] <- sd(R2)
}

# 分位数估计
Q <- quantile((R1 - R_hat)/se_R1, c(1-alpha/2, alpha/2))
names(Q) <- NULL

# R_hat的标准差的bootstrap估计
se_R_hat <- sd(R1)

# t置信区间的bootstrap估计
R_hat - se_R_hat*Q

## -----------------------------------------------------------------------------
library(boot)
# 数据
data <- c(3, 5, 7, 18, 43, 85, 91, 98,100, 130, 230, 487)
# 统计量：均值
boot.mean <- function(x, i) mean(x[i])
# 统计量的bootstrap估计
boot.obj <- boot(data = data, statistic = boot.mean, R = 999)
print(boot.obj)
# 四种区间
ci <- boot.ci(boot.obj, type = c('norm', 'basic', 'perc', 'bca'))
print(ci)

## -----------------------------------------------------------------------------
library(bootstrap)
# 数据集
str(scor)
head(scor)

# theta_hat估计

# 样本协方差阵(为MLE估计) 
sigma <- cor(scor)
# 对样本协方差阵作谱分解
eigen_decompos <- eigen(sigma)
# 求解第一主成分解释的比例
theta_hat <- eigen_decompos[[1]][1]/sum(eigen_decompos[[1]])
theta_hat

# Jacknife方法
n <- nrow(scor)
theta_jack <- numeric(n)
for(i in 1:n){
  sigma_jack <- cor(scor[-i, ])
  theta_jack[i] <- eigen(sigma_jack)[[1]][1]/sum(eigen(sigma_jack)[[1]])
}
# 偏差的Jacknife估计
bias_jack <- (n-1)*(mean(theta_jack) - theta_hat)
bias_jack
# 标准差的Jacknife估计
sd_jack <- sqrt((n-1)*var(theta_jack))
sd_jack

## -----------------------------------------------------------------------------
# 读取并查看数据
library(SA23204171)
data(ironslag)
attach(ironslag)

# leave-two-out cross-validation(n/2 fold CV)
n <- nrow(ironslag)
e1 <- e2 <- e3 <- e4 <- numeric(n)
for(k in 1:(n/2)){
  y <- ironslag$magnetic[-c(2*k-1, 2*k)]
  x <- ironslag$chemical[-c(2*k-1, 2*k)]
  
  J1 <- lm(y~x)
  yhat1 <- J1$coef[1] + J1$coef[2] * ironslag$chemical[c(2*k-1, 2*k)]
  e1[c(2*k-1, 2*k)] <- ironslag$magnetic[c(2*k-1, 2*k)] - yhat1
  
  J2 <- lm(y~x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * ironslag$chemical[c(2*k-1, 2*k)] + 
           J2$coef[3] * ironslag$chemical[c(2*k-1, 2*k)]^2
  e2[c(2*k-1, 2*k)] <- ironslag$magnetic[c(2*k-1, 2*k)] - yhat2
  
  J3 <- lm(log(y)~x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * ironslag$chemical[c(2*k-1, 2*k)]
  yhat3 <- exp(logyhat3)
  e3[c(2*k-1, 2*k)] <- ironslag$magnetic[c(2*k-1, 2*k)] - yhat3
  
  J4 <- lm(log(y)~log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(ironslag$chemical[c(2*k-1, 2*k)])
  yhat4 <- exp(logyhat4)
  e4[c(2*k-1, 2*k)] <- ironslag$magnetic[c(2*k-1, 2*k)] - yhat4
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
# 导入数据
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
str(x); str(y)

n <- length(x); m <- length(y)
R <- 999                            # number of replicates
z <- c(x, y)                        # pooled sample
reps <- numeric(R)                  # storage for replicates

ecdf <- function(sample, x){
  # 经验分布函数
  # sample为对应的样本, 为向量
  # x为经验分布函数的自变量, 可以为向量
  res <- numeric(length(x))
  for(i in 1:length(x)){
    res[i] <- mean(sample <= x[i])
  }
  return(res)
}

W_2 <- function(sample1, sample2){
  # Cramer-von Mises统计量
  n <- length(sample1); m <- length(sample2)
  temp1 <- sum((ecdf(sample1, sample1) - ecdf(sample2, sample1))^2)
  temp2 <- sum((ecdf(sample1, sample2) - ecdf(sample2, sample2))^2)
  return(m*n * (temp1 + temp2)/ (m+n)^2)
}

# Permutation test with Cramer-von Mises statistic
W0 <- W_2(x, y)
for (i in 1:R){
  k <- sample(n+m, size = n, replace = F)
  x1 <- z[k]
  y1 <- z[-k]
  reps[i] <- W_2(x1, y1)
}
p <- mean(c(W0, reps) >= W0); p

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "sunflower"])) # x发生了改变
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
str(x); str(y)

n <- length(x); m <- length(y)
R <- 999                            # number of replicates
z <- c(x, y)                        # pooled sample
reps <- numeric(R)                  # storage for replicates

ecdf <- function(sample, x){
  # 经验分布函数
  # sample为对应的样本, 为向量
  # x为经验分布函数的自变量, 可以为向量
  res <- numeric(length(x))
  for(i in 1:length(x)){
    res[i] <- mean(sample <= x[i])
  }
  return(res)
}

W_2 <- function(sample1, sample2){
  # Cramer-von Mises统计量
  n <- length(sample1); m <- length(sample2)
  temp1 <- sum((ecdf(sample1, sample1) - ecdf(sample2, sample1))^2)
  temp2 <- sum((ecdf(sample1, sample2) - ecdf(sample2, sample2))^2)
  return(m*n * (temp1 + temp2)/ (m+n)^2)
}

# Permutation test with Cramer-von Mises statistic
W0 <- W_2(x, y)
for (i in 1:R){
  k <- sample(n+m, size = n, replace = F)
  x1 <- z[k]
  y1 <- z[-k]
  reps[i] <- W_2(x1, y1)
}
p <- mean(c(W0, reps) >= W0); p

## -----------------------------------------------------------------------------
# count five statistic
count5 <- function(x, y){
  x1 <- x - mean(x)
  y1 <- y - mean(y)
  x_extreme <- sum(x1 > max(y1)) + sum(x1 < min(y1))
  y_extreme <- sum(y1 > max(x1)) + sum(y1 < min(x1))
  return(max(x_extreme, y_extreme))
}

N <- 50       # repeat the test N times
m <- 50; n <- 80
x <- rnorm(m, 0, 1); y <- rnorm(n, 0, 1.5)
D0 <- count5(x, y)
R <- 1000       # number pf replicates
D <- numeric(R) # storage for relicates
p <- numeric(N) # storage for p value in each repeated test

for (i in 1:N){
  for (j in 1:R){
    k <- sample(m+n, size = m, replace = F)
    x1 <- c(x, y)[k]
    y1 <- c(x, y)[-k]
    D[j] <- count5(x1, y1)
  }
  p[i] <- mean(c(D0, D) >= D0)
}
round(p, 4)
mean(p < 0.05)

## -----------------------------------------------------------------------------
intercept <- function(N, b1, b2, b3, f0){
  x1 <- rpois(N, lambda = 1)
  x2 <- rexp(N, rate = 1)
  x3 <- rbinom(N, size = 1, prob = 0.5)
  
  logit <- function(a){
    t <- exp(a + b1*x1 + b2*x2 + b3*x3)
    p <- t / (t + 1)
    return(mean(p) - f0)
  }
  
  res <- uniroot(logit, c(-30, 30))
  return(unlist(res)[1:3])
}

N <- 1e6
b1 <- 0; b2 <- 1; b3 <- -1
f0 <- c(0.1, 0.01, 0.001, 0.0001)

intercept(N, b1, b2, b3, f0[1])
intercept(N, b1, b2, b3, f0[2])
intercept(N, b1, b2, b3, f0[3])
intercept(N, b1, b2, b3, f0[4])

f0 <- seq(0.001, 0.02, 0.001)
a <- numeric(20)
for (i in 1:20){
  a[i] <- intercept(N, b1, b2, b3, f0[i])[1]
}
neg_log_f0 <- -log(f0)
plot(neg_log_f0, a, main = '-log(f0) vs a', col = 'red', pch = 19)

## -----------------------------------------------------------------------------
RW_Metropolis_Laplace <- function(sigma, x0, N){
  # 用random walk metropolis方法生成标准拉普拉斯随机数
  x <- numeric(N) # 生成的随机游走链的长度
  x[1] <- x0 # 用户给定初始值
  u <- runif(N)
  k <- 0 # 标记拒绝的次数
  for (i in 2:N){
    y <- rnorm(1, x[i-1], sigma) # 提议分布采用正态分布，方差由用户指定
    if (u[i] <= (exp(abs(x[i-1]) - abs(y)))){
      x[i] <- y
    } else{
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x = x, rej_num = k, accept_rate = (N-k)/N))
}

x0 <- 10
N <- 2000
sigma <- c(0.05, 0.5, 2, 8)
rw1 <- RW_Metropolis_Laplace(sigma[1], x0, N)
rw2 <- RW_Metropolis_Laplace(sigma[2], x0, N)
rw3 <- RW_Metropolis_Laplace(sigma[3], x0, N)
rw4 <- RW_Metropolis_Laplace(sigma[4], x0, N)

cat('sigma = 0.05: rej_num is ', rw1$rej_num, ', accept_rate is ', rw1$accept_rate, '\nsigma = 0.5 : rej_num is ', rw2$rej_num, ', accept_rate is ', rw2$accept_rate,
'\nsigma = 2   : rej_num is ', rw3$rej_num, ', accept_rate is ', rw3$accept_rate,
'\nsigma = 8   : rej_num is ', rw4$rej_num, ', accept_rate is ', rw4$accept_rate, sep = '')

index <- 1000:1500
plot(index, rw1$x[index], type = 'l', main = 'sigma = 0.05', xlab = '', ylab = '')
plot(index, rw2$x[index], type = 'l', main = 'sigma = 0.5', xlab = '', ylab = '')
plot(index, rw3$x[index], type = 'l', main = 'sigma = 2', xlab = '', ylab = '')
plot(index, rw4$x[index], type = 'l', main = 'sigma = 8', xlab = '', ylab = '')

## -----------------------------------------------------------------------------
N <- 5000            # length of chain
burn <- 2000         # burn-in length
X <- matrix(0, nrow = N, ncol = 2) # the chain, a bivariate sample

ro <- 0.9                          # correlation
mu1 <- 0; mu2 <- 0                 # mean vector
sigma1 <- 1; sigma2 <- 1           # standard deviation
cond_sigma1 <- sqrt(1-ro^2)*sigma1 # conditional variance of x given y
cond_sigma2 <- sqrt(1-ro^2)*sigma2 # conditional variance of y given x

X[1, ] <- c(mu1, mu2) # initialize
for (i in 2:N){
  y <- X[i-1, 2]
  cond_mu1 <- mu1 + ro * (y-mu2) * sigma1/sigma2 # conditional mean of x given y
  X[i, 1] <- rnorm(1, cond_mu1, cond_sigma1)
  x <- X[i, 1]
  cond_mu2 <- mu2 + ro * (x-mu1) * sigma2/sigma1 # conditional mean of y given x
  X[i, 2] <- rnorm(1, cond_mu2, cond_sigma2)
}

# plot the generated sample after discarding a burn-in period of length 2000
b <- burn + 1
x <- X[b:N, ]
plot(x, main = '', cex = 0.25, xlab = 'x', ylab = 'y', xlim = range(x[, 1]), ylim = range(x[, 2]), pch = 19, col = 'red')

# fit a simple linear model
fit <- lm(x[, 2]~x[, 1])
summary(fit)

# check the residuals for normality and constant variance
##### 残差正态性良好 #####
##### 残差无非线性趋势，方差齐性良好 #####
plot(fit)

## -----------------------------------------------------------------------------
# The Rayleigh(sigma) density
f <- function(x, sigma){
  if (any(x < 0)) return(0)
  stopifnot(sigma > 0)
  return((x/sigma^2) * exp(-x^2/(2*sigma^2)))
}

# Metropolis approach to generate a chain of length n
MH_Rayleigh <- function(n, sigma){
  x <- numeric(n)
  x[1] <- rchisq(1, df = 1)
  u <- runif(n)
  
  for (i in 2:n){
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den){
      x[i] <- y # y is accepted
    } else{
      x[i] <- xt # y is rejected
    }
  }
  return(x)
}

# get R_hat
Gelman_Rubin <- function(psi){
  k <- nrow(psi)
  n <- ncol(psi)
  
  psi_means <- rowMeans(psi)
  B <- n * var(psi_means)
  W <- mean(apply(psi, 1, var))
  v_hat <- W*(n-1)/n + B/n
  r_hat <- v_hat/W
  return(r_hat)
}

sigma <- 0.5
k <- 4
b <- 1000
N <- 20000
X <- matrix(0, nrow = k, ncol = N)

# generate chains
for (i in 1:k){
  X[i, ] <- MH_Rayleigh(N, sigma)
}

# compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
rhat <- rep(0, N)
for (i in 1:k){
  psi[i, ] <- psi[i, ]/(1:N)
}
for (j in (b+1):N){
  rhat[j] <- round(Gelman_Rubin(psi[, 1:j]), 3)
}

# plot the sequence of rhat statistics to see when rhat < 1.2
plot(rhat[(b+1):N], type = 'l', xlab = '', ylab = 'Rhat', ylim = c(1, 2))
abline(h = 1.2, lwd = 0.8, col = 'red', lty = 3)

## -----------------------------------------------------------------------------
n <- 10
u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)

# observed data log-likelihood maximization
l_tilta <- function(lambda = 1)  -sum(log(exp(-lambda * u) - exp(-lambda * v)))

library(stats4)
options(warn = -1)
fit <- mle(l_tilta)
fit@coef

# EM algorithm
err <- .Machine$double.eps^0.25 # 收敛所容忍的误差
lambda <- 20                    # EM算法迭代初值
flag <- 1                       # 误差标记
iter <- 0                       # 迭代次数

f <- function(x){
  # 迭代公式
  temp <- sum((u*exp(-x*u) - v*exp(-x*v))/(exp(-x*u) - exp(-x*v)))
  return((n * x/(n + x * temp)))
}

while(flag > err){
  flag <- abs(lambda - f(lambda))
  lambda <- f(lambda)
  iter <- iter + 1
}

cat("EM estimate: ", lambda, "\ntimes of iteration: ", iter, sep = '')

## -----------------------------------------------------------------------------
library(boot)

solve_game <- function(A){
  min_A <- min(A)
  A <- A - min_A
  max_A <- max(A)
  A <- A/max_A
  m <- nrow(A)
  n <- ncol(A)
  it <- n^3
  a <- c(rep(0, m), 1)
  A1 <- -cbind(t(A), rep(-1, n))
  b1 <- rep(0, n)
  A3 <- t(as.matrix(c(rep(1, m), 0)))
  b3 <- 1
  sx <- simplex(a = a, A1 = A1, b1 = b1, A3 = A3, b3 = b3, maxi = T, n.iter = it)
  
  a <- c(rep(0, n), 1)
  A1 <- cbind(A, rep(-1, m))
  b1 <- rep(0, m)
  A3 <- t(as.matrix(c(rep(1, n), 0)))
  b3 <- 1
  sy <- simplex(a = a, A1 = A1, b1 = b1, A3 = A3, b3 = b3, maxi = F, n.iter = it)
  
  soln <- list("A" = A * max_A + min_A, 
               "x" = sx$soln[1:m],
               "y" = sy$soln[1:n],
               "v" = sx$soln[m+1] * max_A + min_A)
  return(soln)
}

# input the payoff matrix
A <- matrix(c(0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
B <- A + 2

# solve game A and game B
s_A <- solve_game(A)
s_B <- solve_game(B)

# solution of game B
round(cbind(s_B$x, s_B$y), 7)

# value of game A and B
value <- c(s_A$v, s_B$v)
names(value) <- c("A", "B")
cat('\n'); print(value)

# one of the extreme points
print(c(0, 0, 25/61, 0, 20/61, 0, 16/61, 0, 0))
# 显然game B的解为game A的一个极值点. 

## -----------------------------------------------------------------------------
# Ex.1
v <- c(1, 2, 3)
dim(v)

# Ex.2
x <- matrix(c(1,1,1,1), 2,2)
is.matrix(x)
is.array(x)

## -----------------------------------------------------------------------------
# Ex.2
d <- data.frame('1'=c(1, 2, 3), '2'=c('a','b','c'), '3'=c(T,T,F), '4'=c(NA,NA,NA))
d
as.matrix(d)

# Ex.3
d_nocol <- d[, -c(1,2,3,4)]
d_nocol
str(d_nocol)
d_norow <- d[-c(1,2,3), ]
d_norow
str(d_norow)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

# apply scale01 to every column of a data frame
d1 <- data.frame('1'=c(2,1,3), '2'=c(4,1,6), '3'=c(5,6,9))
d1
apply(d1, MARGIN = 2, scale01)

# apply scale01 to every numeric column of a data frame
d2 <- data.frame('1'=c(2,1,3), '2'=c(4,1,6), '3'=c('a','b','c'))
d2
cbind(as.data.frame(apply(d2[, as.logical(vapply(d2, is.numeric, 2))], MARGIN = 2, scale01)), d2[, !as.logical(vapply(d2, is.numeric, 2)), drop = F])

## -----------------------------------------------------------------------------
# Compute the standard deviation of every column in a numeric data frame.
d1 <- data.frame('1'=c(2,1,3), '2'=c(4,1,6), '3'=c(5,6,9))
d1
vapply(d1, sd, FUN.VALUE = 1)

# Compute the standard deviation of every numeric column in a mixed data frame. 
d2 <- data.frame('1'=c(2,1,3), '2'=c(4,1,6), '3'=c('a','b','c'))
d2
vapply(d2[, as.logical(vapply(d2, is.numeric, 1))], sd, FUN.VALUE = 1)

## -----------------------------------------------------------------------------
# settings
Nchain <- 500
a <- 1
b <- 2
n <- 20
thin=20

# gibbs in R
gibbsR <- function(N, thin) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0.5
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1, n, y)
      y <- rbeta(1, x+a, n-x+b)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}

