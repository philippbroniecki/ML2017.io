# Code Day 5

## Subset Selection Methods

## the data
rm( list = ls() )

# load foreigners data
load("BSAS_manip_missings.RData")
df <- data2
rm(data2)


# check for missing values
apply(df, 2, function(x) table(is.na(x))["TRUE"] )


# we drop missings in IMMBRIT
df <- df[ !is.na(df$IMMBRIT), ]
df <- na.omit(df)


# we declare the factor variables
df$urban <- factor(df$urban, labels = c("rural", "more rural", "more urban", "urban"))
df$RSex <- factor(df$RSex, labels = c("Male", "Female"))
df$health.good <- factor(df$health.good, labels = c("bad", "fair", "fairly good", "good") )


## best subset
library(leaps)

# drop the binary response coded 1 if IMMBRIT > 10.7 
df$over.estimate <- NULL
df$Cons <- NULL

# run best subset selection
regfit.full <- regsubsets(IMMBRIT ~ ., data = df)
summary(regfit.full)

# increase the max number of variables (defaults to 8)
regfit.full <- regsubsets(IMMBRIT ~ ., data = df, nvmax = 16)
reg.summary <- summary(regfit.full)

# check components of summary object
names(reg.summary)

# R^2's for best model for the number of variables included
reg.summary$rsq

# visualize
par( mfrow =  c(2,2) ) # 2 row, 2 columns in plot window
plot(reg.summary$rss, xlab = "Number of Variables", ylab = "Residual Sum of Squares", type = "l")
plot(reg.summary$adjr2, xlab = "Number of Variables", ylab = "Adjusted R^2", type = "l")

# find the peak of adj. R^2
adjr2.max <- which.max( reg.summary$adjr2 )
points(adjr2.max, reg.summary$adjr2[adjr2.max], col = "red", pch = 20, cex = 2)


# cp
plot(reg.summary$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
cp.min <- which.min(reg.summary$cp)
points(cp.min, reg.summary$cp[cp.min], col = "red", cex = 2, pch = 20)


# bic
bic.min <- which.min(reg.summary$bic)
plot(reg.summary$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")
points(bic.min, reg.summary$bic[bic.min], col = "red", cex = 2, pch = 20)


# plot model comparison based on R^2
par( mfrow = c(1,1) )
plot(regfit.full, scale = "r2")


# plot model comparison based on adjusted R^2
par( mfrow = c(1,1) )
plot(regfit.full, scale = "adjr2")


# plot model comparison based on adjusted CP
par( mfrow = c(1,1) )
plot(regfit.full, scale = "Cp")


# plot model comparison based on adjusted BIC
par( mfrow = c(1,1) )
plot(regfit.full, scale = "bic")


# model coeffiecients of model with min. BIC
coef(regfit.full, id = bic.min)
coef(regfit.full, id = 5)


#### Forward and Backward Stepwise Selection

# run forward selection
regfit.fwd <- regsubsets(IMMBRIT ~ ., data = df, nvmax = 16, method = "forward")
summary(regfit.fwd)

# run backward selection
regfit.bwd <- regsubsets(IMMBRIT ~ ., data = df, nvmax = 16, method = "backward")
summary(regfit.bwd)

# model coefficients of best 7-variable models
coef(regfit.full, 7)
coef(regfit.fwd, 7)
coef(regfit.bwd, 7)


### Choosing Among Models Using the Validation Set Approach and Cross-Validation

# initialize random number generator
set.seed(1)

# sample true or false for each observation
train <- sample( c(TRUE, FALSE), size = nrow(df), replace = TRUE )
# the complement
test <- (!train)


# best subset selection on training set only
regfit.best <- regsubsets(IMMBRIT ~ ., data = df[train, ], nvmax = 16)

# test data
test.mat <- model.matrix(IMMBRIT ~., data = df[test, ])

# validation error for each model
val.errors <- NA
for (i in 1:16 ){
  coefi <- coef(regfit.best, id = i)
  y_hat <- test.mat[, names(coefi)] %*% coefi
  val.errors[i] <- mean(  (df$IMMBRIT[test] - y_hat)^2   )
}

# examine errors
val.errors

# which model has smallest error
min.val.errors <- which.min(val.errors)
# coefficients of that model
coef( regfit.best, min.val.errors )


# precict function for repeatedly choosing model with lowest test error
predict.regsubsets <- function( object, newdata, id, ... ){
  m.formula <- as.formula( object$call[[2]] )
  mat <- model.matrix( m.formula, newdata )
  coefi <- coef( object, id = id )
  xvars <- names( coefi )
  mat[ , xvars ] %*% coefi
}


# best subset on full data set
regfit.best <- regsubsets( IMMBRIT ~ ., data = df, nvmax = 16 )

# examine coefficients of the model that had the lowest validation error
coef( regfit.best, min.val.errors )



## k-fold cross-validation

# number of folds
k <- 10
set.seed(1)
# fold assignment for each observation
folds <- sample(1:k, nrow(df), replace = TRUE)
table(folds)

# container for cross-validation errors
cv.errors <- matrix(NA, nrow = k, ncol = 16, dimnames = list(NULL, paste(1:16)))
cv.errors

# loop over folds
for (a in 1:k){
  
  # best subset selection on training data
  best.fit <- regsubsets(IMMBRIT ~ ., data = df[ folds != a, ], nvmax = 16)
  
  # loop over the 16 subsets
  for (b in 1:16){

    # predict response for test set for current subset
    pred <- predict(object = best.fit, newdata = df[ folds == a ,], id = b )
    
    # MSE into container; rows are folds; columns are subsets
    cv.errors[a, b] <- mean( (df$IMMBRIT[folds==a] - pred)^2 )
    
  } # end of loop over the 20 subsets
} # end of loop over folds
cv.errors

# average cross-validation errors over the folds
mean.cv.errors <- apply(cv.errors, 2, mean)
mean.cv.errors

# visualize
par( mfrow = c(1,1) )
plot( mean.cv.errors, type = "b" )


# run regsubsets on full data set
reg.best <- regsubsets(IMMBRIT ~ ., data = df, nvmax = 16)
# coefficients of subset which minimized test error
coef(reg.best, which.min(mean.cv.errors))

