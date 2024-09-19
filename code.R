setwd("/projectnb/econdept/estamm/EC709_PS1/")

## Load dependencies, install if not already.
packages <-
  c('splines','locfit','dplyr',
    'parallel',
    'foreach',
    'doParallel','hdm','glmnet','mvtnorm',
    'hdm')

for (pkg in packages) {
  if (require(pkg, character.only = TRUE) == FALSE) {
    print(paste0("Trying to install ", pkg))
    install.packages(pkg)
    if (require(pkg, character.only = TRUE)) {
      print(paste0(pkg, " installed and loaded"))
    } else{
      stop(paste0("could not install ", pkg))
    }
  }
}

# Question 1.3----
set.seed(1004)
x <- 0.5
n <- 1000
h <- 1:4*0.05
X <- runif(n)
e <- rnorm(n)
Y <- exp(X)*(1+e)

results <- data.frame(h = h, kern_bias = rep(NA,4), 
                      kern_var = rep(NA,4), local_bias = rep(NA, 4),
                      local_var = rep(NA, 4))

for(i in h){
  kernreg <- locpoly(x = X, y = Y, bandwidth = i,
                     range.x = c(0,1), gridsize = 500,
                     degree = 0)
  locreg <- locpoly(x = X, y = Y, bandwidth = i,
                    range.x = c(0,1), gridsize = 500,
                    degree = 1)
  results$kern_bias[results$h == i] <- mean(kernreg$y - Y)
  results$local_bias[results$h == i] <- mean(locreg$y - Y)
  results$kern_var[results$h == i] <- var(kernreg$y)
  results$local_var[results$h == i] <- var(locreg$y)
  
}

stargazer(results, type = "latex", out = "q13.tex", summary = FALSE)
# Q2
df <- read_dta("yn.dta")
df$lprice <- log(df$price)
df$lgas <- log(df$gas)
df$lincome <- log(df$income)


crossvalidate <- function(h,X,Y,type){
  n <- nrow(X)
  nx <- ncol(X)
  
  loss <- numeric(n)
  for(i in 1:n){
    X_tr <- X[-i,]
    Y_tr <- Y[-i,]
    X_test <- X[i,]
    Y_test <- Y[i,]
    train <- data.frame(Y_tr,X_tr)
    test <- data.frame(Y_test,X_test)
    colnames(train) <-  c("Y",paste0("X",1:nx))
    colnames(test) <-  c("Y",paste0("X",1:nx))
    
    if(type == "kern"){
      
      ff <- as.formula(paste("Y ~ lp(", paste0("X",1:nx, collapse = ", "),", degree = 0, h = ",h,")"))
      mod2 <- locfit(data = train, ff)
      X_test_final <- test %>% select(-Y)
      
      Y_pred <- predict(mod, newdata =X_test_final)  
      
    }
    else if(type == "loc"){
      ff <- as.formula(paste("Y ~ lp(", paste0("X",1:nx, collapse = ", "),", degree = 1, h = ",h,")"))
      
      mod <- locfit(data = train, ff)
      X_test_final <- test %>% select(-Y)
      
      Y_pred <- predict(mod,newdata =X_test_final)
      
    }
    else if(type == "power"){
      # This shit is problematic so using a workaround stolen from here: https://stackoverflow.com/questions/34109923/r-predicting-from-multivariate-polynomial-models
      
      ff <- as.formula(paste("Y ~ poly(", paste0("X",1:nx, collapse = ", "),", degree = ",h,")" ))
      mod <- lm(data = train, ff)
      X_test_final <- rbind(select(test,-Y),c(0,0))
      Y_pred <- predict(mod,newdata = X_test_final)[1]
      
    }
    else if(type == "bspline"){
      
      # This shit is problematic so using a workaround stolen from here: https://stackoverflow.com/questions/34109923/r-predicting-from-multivariate-polynomial-models
      
      ff <- as.formula(paste("Y ~ ", paste0("bs(X",1:nx,", knots = ",h,")", collapse = " + ")))
      mod <- lm(data = train, ff)
      X_test_final <- rbind(select(test,-Y),c(0,0))
      Y_pred <- predict(mod,newdata = X_test_final)[1]
      
    }
    loss[i] <- (Y_test[[1]][1] - Y_pred)^2
    
  }
  return(mean(loss))
}

h_kern <- 1:8*0.1
h_loc <- 1:8*0.1
h_power <- 1:7
h_bspline <- 1:7

X <- df %>% select(lprice)
Y <- df %>% select(lgas)

ncore <- detectCores()
cl <- makeCluster(ncore - 1)
registerDoParallel(cl)

locloss_2_1<- foreach(a=h_loc, .packages = c('splines','locfit','dplyr'), .verbose = TRUE, .combine = rbind) %dopar% {
  crossvalidate(a,X,Y,'loc')
}

powerloss_2_1<- foreach(a=h_power, .packages = c('splines','locfit','dplyr'), .verbose = TRUE, .combine = rbind) %dopar% {
  crossvalidate(a,X,Y,'power')
}

bsplineloss_2_1<- foreach(a=h_bspline, .packages = c('splines','locfit','dplyr'), .verbose = TRUE, .combine = rbind) %dopar% {
  crossvalidate(a,X,Y,'bspline')
}


locmin_2_1 <- which.min(locloss_2_1)
powermin_2_1 <- which.min(powerloss_2_1)
bsplinemin_2_1 <- which.min(bsplineloss_2_1)

loc_2_1 <- locfit(data = df, lgas ~ lp(lprice, degree = 1, bandwidth = h_loc[locmin_2_1]))
power_2_1 <- lm(data = df,lgas ~ poly(lprice, degree = h_power[powermin_2_1]))
bspline_2_1 <- lm(data = df,lgas ~ bs(lprice, knots = h_bspline[bsplinemin_2_1]))

# Create a sequence of lprice values for prediction
lprice_pred <- seq(min(df$lprice), max(df$lprice), length.out = 100)

# Predict lgas values
lgas_pred_loc <- predict(loc_2_1, newdata = data.frame(lprice = lprice_pred))
lgas_pred_power <- predict(power_2_1, newdata = data.frame(lprice = lprice_pred))
lgas_pred_bspline <- predict(bspline_2_1, newdata = data.frame(lprice = lprice_pred))

# Create the plot
png('Q2_1.png')
plot(df$lprice, df$lgas, main="lgas ~ lp(lprice)", 
     xlab="lprice", ylab="lgas", pch=20, col="gray")
lines(lprice_pred, lgas_pred_loc, col="red", lwd=2)
lines(lprice_pred, lgas_pred_power, col="blue", lwd=2)
lines(lprice_pred, lgas_pred_bspline, col="green", lwd=2)
legend("bottomright", legend=c("Local linear",'Power Series','BSpline'), 
       col=c("red", "blue",'green'), lty=c(1,2), lwd=c(2,1))
dev.off()

X_new <- log(0.57)

Y_pred_loc <- predict(loc_2_1, se.fit = TRUE, data.frame(lprice = X_new))
Y_pred_power <- predict(power_2_1, se.fit = TRUE,data.frame(lprice = X_new))
Y_pred_bspline <- predict(bspline_2_1, se.fit = TRUE,data.frame(lprice = X_new))




X <- df %>% select(lprice,lincome)
Y <- df %>% select(lgas)

locloss_2_3<- foreach(a=h_loc, .packages = c('splines','locfit','dplyr'), .verbose = TRUE, .combine = rbind) %dopar% {
  crossvalidate(a,X,Y,'loc')
}

powerloss_2_3<- foreach(a=h_power, .packages = c('splines','locfit','dplyr'), .verbose = TRUE, .combine = rbind) %dopar% {
  crossvalidate(a,X,Y,'power')
}

bsplineloss_2_3<- foreach(a=h_bspline, .packages = c('splines','locfit','dplyr'), .verbose = TRUE, .combine = rbind) %dopar% {
  crossvalidate(a,X,Y,'bspline')
}


stopCluster(cl)

powermin_2_3 <- which.min(powerloss_2_3)
bsplinemin_2_3 <- which.min(bsplineloss_2_3)
locmin_2_3 <- which.min(locloss_2_3)

loc_2_3 <- locfit(data = df, lgas ~ lp(lprice, lincome, degree = 1, bandwidth = h_loc[locmin_2_3]))
power_2_3 <- lm(data = df,lgas ~ poly(lprice, lincome, degree = h_power[powermin_2_3]))
bspline_2_3 <- lm(data = df,lgas ~ bs(lprice, knots = h_bspline[bsplinemin_2_3]) +  bs(lincome, knots = h_bspline[bsplinemin_2_3]))

#2.4
lasso <- rlasso(data = df ,lgas ~ lprice + age + driver + hhsize + fueltype + urban + prov + year , post = FALSE)  # use lasso, not-Post-lasso

lasso_twoway <- rlasso(data = df ,lgas ~ (lprice + age + driver + hhsize + fueltype + urban + prov + year)^2 , post = FALSE)  # use lasso, not-Post-lasso


# Q3
set.seed(859)
n <- 100
K <- 500
sims <- 500
b_0 <- c(1,1,0.5,2/3,1/4,1/5,rep(0,K-6))
n_folds <- 10
alpha <- 1

pred_error <- function(bhat,X,Xb_0){
  Xbhat <- X %*% bhat
  mean((Xbhat - Xb_0)^2)
}

lasso_pe <- numeric(sims)
post_lasso_pe <- numeric(sims)
lassoCV_pe <- numeric(sims)
post_lassoCV_pe <- numeric(sims)
oracle_pe <- numeric(sims)
Sigma <- outer(1:K,1:K, function(i,j) 0.5^abs(i - j))

for(i in 1:sims){
  lasso_betas <- numeric(K)
  post_lasso_betas <- numeric(K)
  lassoCV_betas <- numeric(K)
  post_lassoCV_betas <- numeric(K)
  oracle_betas <- numeric(K)
  
  e <- rnorm(n)
  
  X <- rmvnorm(n, sigma = Sigma)
  Xb_0 <- X %*% b_0
  Y <- Xb_0 + e
  
  lasso <- glmnet(X,Y, alpha = alpha, intercept = FALSE, nlambda = 1)
  lasso_coefs <- coef(lasso)
  lasso_betas <- lasso_coefs[-1]
  lasso_pe[i] <- pred_error(lasso_betas,X,Xb_0)
  
  post_lasso <- lm(Y ~ 0 + X[,lasso_coefs@i[lasso_coefs@i != 0]])
  post_lasso_betas[lasso_coefs@i[lasso_coefs@i != 0]] <- coef(post_lasso)
  post_lasso_pe[i] <- pred_error(post_lasso_betas,X,Xb_0)
  
  lassoCV <- cv.glmnet(X,Y,alpha = alpha, nfolds = n_folds, intercept = FALSE)
  lassoCV_coefs<- coef(lassoCV, s="lambda.min")
  lassoCV_betas <- lassoCV_coefs[-1]
  lassoCV_pe[i] <- pred_error(lassoCV_betas,X,Xb_0)
  
  
  post_lassoCV <- lm(Y ~ 0 + X[,lassoCV_coefs@i[lassoCV_coefs@i != 0]])
  post_lassoCV_betas[lassoCV_coefs@i[lassoCV_coefs@i != 0]] <- coef(post_lassoCV)
  post_lassoCV_pe[i] <- pred_error(post_lassoCV_betas,X,Xb_0)
  
  oracle <- lm(Y ~ 0 + X[,1:6])
  oracle_betas <- c(coef(oracle), rep(0,K-6))
  oracle_pe[i] <- pred_error(oracle_betas,X,Xb_0)
}

results <- data.frame(lasso_pe,post_lasso_pe,lassoCV_pe,post_lassoCV_pe, oracle_pe)
stargazer(results, type = "latex", out = "q3.tex")
