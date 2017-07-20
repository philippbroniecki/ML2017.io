# Lab Day 9

rm(list = ls())

# data on violent crimes by US state
library(MASS)
?USArrests
head(USArrests)

# create a state variable
states <- row.names(USArrests)
states

# variable names in the data set
names(USArrests)

# summary stats mean and variance
apply(USArrests, 2, mean)
apply(USArrests, 2, var)

#_________________________________________________________________
# pca
#_________________________________________________________________

# perform pca on data
pr.out <- prcomp(USArrests, scale = TRUE)

# check model object
names(pr.out)

pr.out$center
pr.out$scale

# "translation from x to z..."
pr.out$rotation

# dimensions of facorized output
dim(pr.out$x)
# the factor scores for each observation
head(pr.out$x) # this what we could use as IVs

# which component picks up most of the variance on the variables
biplot(pr.out, scale = 0)

# turning things around
pr.out$rotation <- -pr.out$rotation
pr.out$x <- -pr.out$x
biplot(pr.out, scale = 0)

# standard deviation on new scales
pr.out$sdev

# variance
pr.var <- pr.out$sdev^2
pr.var

# amount of variance explained by each component
pve <- pr.var / sum(pr.var)
pve

# Finding the ellbow
plot(pve, xlab = "Principal Component", 
ylab = "Proportion of Variance Explained ", 
ylim = c(0, 1), type = "b")

# new command - cumulative sums
a <- c(1, 2, 8, -3)
cumsum(a)

# how much variance do we explain with each additional pc? What if # pc's = # vars?
plot(cumsum(pve), xlab = "Principal Component ",
ylab = " Cumulative Proportion of Variance Explained ", 
ylim = c(0, 1), type = "b")




#________________________________________________________________________________________
## Clustering
#________________________________________________________________________________________

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### K-Means Clustering
set.seed(2)

# fake data; 2 columns of 50 obs from standard normal
x <- matrix(rnorm(50 * 2), ncol = 2)
# adding a systematic offset to first half of col 1 and 2
x[1:25, 1] <- x[1:25, 1] + 3
x[1:25, 2] <- x[1:25, 2] - 4

# visualize the data
plot(x, pch = 16, bty = "n")

# run k-means clustering
km.out <- kmeans(x, 2, nstart = 20)
# for each observation what cluster has it been assigned to?
km.out$cluster

# how where the observations clustered
plot(x, col = (km.out$cluster + 1), # plus b/c the first color is black 
main = "K-Means Clustering Results with K=2", 
xlab = "", ylab = "", pch = 20, cex = 2, bty = "n")

# k=3
set.seed(4)
km.out <- kmeans(x, 3, nstart = 20)
km.out

# visualize cluster assignment again
plot(x, col = (km.out$cluster + 1), 
main = "K-Means Clustering Results with K=3", 
xlab = "", ylab = "", pch = 20, cex = 2)

# set nstart large enough
set.seed(3)
km.out <- kmeans(x, 3, nstart = 1)
km.out$tot.withinss
km.out <- kmeans(x, 3, nstart = 20)
km.out$tot.withinss

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Hierarchical Clustering

# how to estimate which is the next closest observation
hc.complete <- hclust(dist(x), method = "complete")
hc.average <- hclust(dist(x), method = "average")
hc.single <- hclust(dist(x), method = "single")

# complete (maximum distance)
par(mfrow = c(1, 3))
plot(hc.complete, main = "Complete Linkage", 
xlab = "", sub = "", cex = 0.9)

# average
plot(hc.average, main = "Average Linkage", 
xlab = "", sub = "", cex = 0.9)

# mimimum distance
plot(hc.single, main = "Single Linkage", 
xlab = "", sub = "", cex = 0.9)

# cluster assignment for 2 clusters
cutree(hc.complete, 2)
cutree(hc.average, 2)
cutree(hc.single, 2)

# cluster assignment for 4 clusters
cutree(hc.single, 4)

# scaling to get variable on the same scale
xsc <- scale(x)
par( mfrow = c(1,1) )
plot(hclust(dist(xsc), method = "complete"), 
main = "Hierarchical Clustering with Scaled Features ")

# use different similarity measure
x <- matrix(rnorm(30 * 3), ncol = 3) # new fake data
dd <- as.dist(1 - cor(t(x) ))
plot(hclust(dd, method = "complete"), 
main = "Complete Linkage with Correlation-Based Distance", 
xlab = "", sub = "")



## Example from James et al
# gene expression data
library(ISLR)
nci.labs <- NCI60$labs
nci.data <- NCI60$data
?NCI60
dim(nci.data)

# we check the data
nci.data[1:5, 1:5]
head(nci.labs)


### 9.3.1 PCA on the NCI60 Data
pr.out <- prcomp(nci.data, scale = TRUE)

# make a nice plot - need a function
Cols <- function(vec) {
  cols <- rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

# plot
par(mfrow = c(1, 2))
plot(pr.out$x[, 1:2], col = Cols(nci.labs), 
     pch = 19, xlab = "Z1", ylab = "Z2")

plot(pr.out$x[, c(1, 3)], col = Cols(nci.labs), 
     pch = 19, xlab = "Z1", ylab = "Z3")

summary(pr.out)
plot(pr.out)


pve <- 100 * pr.out$sdev^2 / sum(pr.out$sdev^2)
par(mfrow = c(1, 2))
plot(pve, type = "o", ylab = "PVE", 
     xlab = "Principal Component", col = " blue ")
plot(cumsum(pve), type = "o", ylab = "Cumulative PVE", 
     xlab = "Principal Component ", col = " brown3 ")



### 9.3.2 Clustering the Observations of the NCI60 Data
sd.data <- scale(nci.data)

par(mfrow = c(1, 3))
data.dist <- dist(sd.data)
plot(hclust(data.dist), labels = nci.labs, 
     main = "Complete Linkage", xlab = "", sub = "", ylab = "")
plot(hclust(data.dist, method = "average"), 
     labels = nci.labs, main = "Average Linkage", 
     xlab = "", sub = "", ylab = "")
plot(hclust(data.dist, method = "single"), 
     labels = nci.labs, main = "Single Linkage", 
     xlab = "", sub = "", ylab = "")


# cluster object
hc.out <- hclust(data.dist)  
hc.clusters <- cutree(hc.out, 4)
table(hc.clusters, nci.labs)

# plot dendrogram
par(mfrow = c(1, 1))
plot(hc.out, labels = nci.labs)
abline(h = 139, col = "red")

# K-means
set.seed(2)
km.out <- kmeans(sd.data, 4, nstart = 20)
km.clusters <- km.out$cluster
# do we get the same clusters?
table(km.clusters, hc.clusters)

