#####################################
# BST 227 Final Project             #
# 2016 Help Guide                   #
# Author: caleblareau@g.harvard.edu #
#####################################

# THIS VERSION EDITED BY ANDY SHI <andyshi@g.harvard.edu>

library(stats)
library(ggplot2)
library(qqman)
library("data.table")
library("dplyr")
library("magrittr")
library("irlba")

#Read data; consider using the `readr` package to make this faster
gdat <- data.matrix(fread("genotype.txt"))
pdat <- fread("phenotype.txt", header = TRUE)
legend <- fread("legend.txt", header = TRUE)

#################
#Quality Control#
#################

#Count individuals and variants with missingness ratios; consider filtering them
snpMissingRate <- apply(gdat, 2, function(x) mean(is.na(x)))
personMissingRate <- apply(gdat, 1, function(x) mean(is.na(x)))

#Remove badly missing data (e.g. at 5%)
filtering_lvl <- 0.05
gdat2 <- gdat[personMissingRate <= filtering_lvl,
              snpMissingRate <= filtering_lvl]
# similarly filter the other data
pdat2 <- pdat[personMissingRate <= filtering_lvl,]
legend2 <- legend[snpMissingRate <= filtering_lvl,]

#After quality control, numbers remaining
nPeopleRemaining <- dim(gdat2)[1]
nSNPsRemaining <- dim(gdat2)[2]

#Count Number of Individuals Belonging to Each Group after QC
pdat2 %>% group_by(Y, gender, ancestry) %>%
    summarize(total = n(), percentage = n() / nPeopleRemaining)

#feCases <-nrow(subset(pdat,gender=="F" & Y==1 & ancestry == "European"))

# further QC: Hardy Weinberg Equilibrium

#############################
#PrincipalComponent Analysis#
#############################

# Consider taking a random subset of 'r' variants to speed up PCA
r <- 1000
gRand <- gdat[,sample(ncol(gdat), r) ]
G <- as.matrix(gRand)

# Replace missing data for PC to work
# need to fix r or gRand
for (j in 1:r) {
    temp <- gRand[, j]
    # can use mean or median
    G[, j][is.na(temp)] <- mean(temp, na.rm = T)
}

# recommend to center and scale the data first

# Calculate principal components using stats package
# Multiple functions that allow us to do this in R


#####################
#Logistic Regression#
#####################

for (j in 1:nSNPsRemaining) {
    if (j %in% floor(quantile(seq(nSNPsRemaining),seq(0,1,0.1)))) {
        print(paste("Progress:", names(which(j==floor(quantile(seq(nSNPsRemaining),seq(0,1,0.1)))))), quote=F)
    }
    noadj.mod <- summary(glm(pdat$Y~gdat[,j]+pdat$gender,family=binomial()))$coeff

}


#######
#Plots#
#######

manhattan(result,chr="chr",bp="position",p="noadj.p", main="No adjustment", suggestiveline = -log10(.05/nSNPsRemaining))
qq(result$noadj.p, main="No adjustment")
