#####################################
# BST 227 Final Project             #
# Team Morton                       #
#####################################

library(stats)
library(ggplot2)
library(qqman)
library("data.table")
library("dplyr")
library("magrittr")
library("irlba")
library("foreach")
library("doParallel")
library("ggplot2")
library(aod)

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

#Count Number of Individuals Belonging to Each Group before QC:
pdat %>% group_by(Y, gender, ancestry) %>%
    summarize(total = n(),percentage = n() / length(pdat$Y))
#Count Number of Individuals Belonging to Each Group after QC:
pdat2 %>% group_by(Y, gender, ancestry) %>%
    summarize(total = n(), percentage = n() / nPeopleRemaining)

# do something
#feCases <-nrow(subset(pdat,gender=="F" & Y==1 & ancestry == "European"))

# further QC: Hardy Weinberg Equilibrium
# x: a vector of genotypes (0, 1, or 2)
hwe_test <- function(x) {
    x_notna <- x[!is.na(x)]
    stopifnot(is.integer(x_notna))
    n <- length(x_notna)
    pbar <- sum(x_notna) / (2 * n)

    obs <- vapply(0:2, function(y) { sum(x_notna == y)}, as.integer(1))
    expected <- c(n * (1 - pbar)^2, 2 * n * pbar * (1 - pbar),
                  n * (pbar^2))
    test_stat <- sum( ((obs - expected)^2) / expected)
    pval <- pchisq(test_stat, df = 1, lower.tail = FALSE)
    return(pval)
}

hwe_res <- apply(gdat2, 2, hwe_test)
hist(hwe_res, main = "Histogram of p-values",
     xlab = "pvalue for test for HWE")

# hardy-weinberg by ethnicity
hwe_african <- apply(gdat2[pdat2$ancestry == "African", ], 2, hwe_test)
hwe_euro <- apply(gdat2[pdat2$ancestry == "European", ], 2, hwe_test)

hist(hwe_african, main = "Histogram of p-values for Africans",
     xlab = "pvalue for test for HWE")
hist(hwe_euro, main = "Histogram of p-values for Europeans",
     xlab = "pvalue for test for HWE")

#############################
#PrincipalComponent Analysis#
#############################

# Consider taking a random subset of 'r' variants to speed up PCA
G <- as.matrix(gdat2)

# Replace missing data for PC to work
# need to fix r or gRand
for (j in 1:ncol(G)) {
    temp <- G[, j]
    # can use mean or median
    G[, j][is.na(temp)] <- mean(temp, na.rm = T)
    G[, j] <- (G[, j] - mean(G[, j])) / sd(G[, j])
}

# Calculate principal components using stats package
# recommend to center and scale the data first
# Multiple functions that allow us to do this in R

N_PCS <- 10
#col_means <- colMeans(G)
#col_sds <- apply(G, 2, sd)
start_time <- Sys.time()
#svd_res <- irlba(G, nv = N_PCS, center = col_means,
#                 scale = col_sds)
svd_res <- irlba(G, nv = N_PCS)
elapsed <- Sys.time() - start_time
print(elapsed)

g_pca <- G %*% svd_res$v
colnames(g_pca) <- paste0("PC", 1:N_PCS)
g_pca_df <- as.data.frame(g_pca)

ancestry_gender <- interaction(pdat2$ancestry, pdat2$gender)

ggplot(data = g_pca_df, aes(x = PC1, y = PC2, color = pdat2$ancestry)) +
    geom_point() + labs(color = "Ancestry") +
    ggtitle("Principal Components Plot, PC1 and PC2")

ggplot(data = g_pca_df, aes(x = PC2, y = PC3, color = pdat2$ancestry)) +
    geom_point() + labs(color = "Ancestry") +
    ggtitle("Principal Components Plot, PC2 and PC3")

plot(svd_res$d, type = "b", main = "PCA Variance Explained",
     xlab = "Component", ylab = "Singular Value")


#####################
#Logistic Regression#
#####################

## pdat3 takes pdat2 and adds the PC values ##
pdat3 = pdat2
pdat3$PC1 = g_pca_df$PC1
pdat3$PC2 = g_pca_df$PC2
pdat3$PC3 = g_pca_df$PC3
pdat3$PC4 = g_pca_df$PC4
pdat3$PC5 = g_pca_df$PC5
pdat3$PC6 = g_pca_df$PC6
pdat3$PC7 = g_pca_df$PC7
pdat3$PC8 = g_pca_df$PC8
pdat3$PC9 = g_pca_df$PC9
pdat3$PC10 = g_pca_df$PC10


cl <- makeCluster(4)
registerDoParallel(cl)

output.gender <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
    mod <- glm(pdat2$Y ~ gdat2[, j] + pdat2$gender,
                family = binomial())
    pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
    beta <- coef(mod)["gdat2[, j]"]
    ret <- c(j, beta, pval)
    names(ret) <- c("SNP_num", "beta", "pval")
    ret
}
stopCluster(cl)
df.gender <- data.frame("chr"=legend[output.gender[,"SNP_num"],"chr"],
                        "bp"=legend[output.gender[,"SNP_num"],"position"],
                        "p"=output.gender[,"pval"])


cl <- makeCluster(4)
registerDoParallel(cl)

output.noadj <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  mod <- glm(pdat2$Y ~ gdat2[, j],
             family = binomial())
  pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
  beta <- coef(mod)["gdat2[, j]"]
  ret <- c(j, beta, pval)
  names(ret) <- c("SNP_num", "beta", "pval")
  ret
}
stopCluster(cl)
df.noadj <- data.frame("chr"=legend[output.noadj[,"SNP_num"],"chr"],
                        "bp"=legend[output.noadj[,"SNP_num"],"position"],
                        "p"=output.noadj[,"pval"])

cl <- makeCluster(4)
registerDoParallel(cl)

output.anc <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  mod <- glm(pdat2$Y ~ gdat2[, j] + pdat2$ancestry,
             family = binomial())
  pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
  beta <- coef(mod)["gdat2[, j]"]
  ret <- c(j, beta, pval)
  names(ret) <- c("SNP_num", "beta", "pval")
  ret
}
stopCluster(cl)
df.anc <- data.frame("chr"=legend[output.anc[,"SNP_num"],"chr"],
                       "bp"=legend[output.anc[,"SNP_num"],"position"],
                       "p"=output.anc[,"pval"])

cl <- makeCluster(4)
registerDoParallel(cl)

output.genanc <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  mod <- glm(pdat2$Y ~ gdat2[, j] + pdat2$gender + pdat2$ancestry,
             family = binomial())
  pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
  beta <- coef(mod)["gdat2[, j]"]
  ret <- c(j, beta, pval)
  names(ret) <- c("SNP_num", "beta", "pval")
  ret
}
stopCluster(cl)
df.genanc <- data.frame("chr"=legend[output.genanc[,"SNP_num"],"chr"],
                     "bp"=legend[output.genanc[,"SNP_num"],"position"],
                     "p"=output.genanc[,"pval"])

## Ancestry+Gender, by factor ##
cl <- makeCluster(4)
registerDoParallel(cl)

output.anc.fac <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  use <- ifelse(is.na(gdat2[,j]),FALSE,TRUE)
  numterms <- sum(c(sum(gdat2[use,j]==0),sum(gdat2[use,j]==1),
                    sum(gdat2[use,j]==2))>0)
  if(numterms >= 2) {
    mod <- glm(pdat2$Y[use] ~ as.factor(gdat2[use, j]) + 
                 pdat2$ancestry[use] + pdat2$gender[use],
               family = binomial())
    wald <- tryCatch(aod::wald.test(Sigma=vcov(mod),b=coef(mod),Terms=2:numterms),
                     error=function(e) aod::wald.test(Sigma=1,
                                                      b=0,
                                                      Terms=1))
    #wald <- try(aod::wald.test(Sigma=vcov(mod),b=coef(mod),Terms=2:numterms))
    #if (inherits(wald, "try-error")) {
    #  pval <- 1
    #}
    pval <- wald$result$chi2["P"]
    ret <- c(j, pval)
  }
  else {
    ret <- c(j,1)
  }
names(ret) <- c("SNP_num", "pval")
ret
}
stopCluster(cl)
df.anc.fac <- data.frame("chr"=legend[output.anc.fac[,"SNP_num"],"chr"],
                     "bp"=legend[output.anc.fac[,"SNP_num"],"position"],
                     "p"=output.anc.fac[,"pval"])


## PC1-3 ##
cl <- makeCluster(4)
registerDoParallel(cl)

output.pc13 <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  mod <- glm(pdat3$Y ~ gdat2[, j] + pdat3$PC1 + pdat3$PC2 + pdat3$PC3,
             family = binomial())
  pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
  beta <- coef(mod)["gdat2[, j]"]
  ret <- c(j, beta, pval)
  names(ret) <- c("SNP_num", "beta", "pval")
  ret
}
stopCluster(cl)
df.pc13 <- data.frame("chr"=legend[output.pc13[,"SNP_num"],"chr"],
                        "bp"=legend[output.pc13[,"SNP_num"],"position"],
                        "p"=output.pc13[,"pval"])

## PC1-5 ##
cl <- makeCluster(4)
registerDoParallel(cl)

output.pc15 <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  mod <- glm(pdat3$Y ~ gdat2[, j] + pdat3$PC1 + pdat3$PC2 + pdat3$PC3 + pdat3$PC4 + pdat3$PC5,
             family = binomial())
  pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
  beta <- coef(mod)["gdat2[, j]"]
  ret <- c(j, beta, pval)
  names(ret) <- c("SNP_num", "beta", "pval")
  ret
}
stopCluster(cl)
df.pc15 <- data.frame("chr"=legend[output.pc15[,"SNP_num"],"chr"],
                      "bp"=legend[output.pc15[,"SNP_num"],"position"],
                      "p"=output.pc15[,"pval"])

## Gender+PC1-3 ##
cl <- makeCluster(4)
registerDoParallel(cl)

output.genpc13 <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  mod <- glm(pdat3$Y ~ gdat2[, j] + pdat3$gender + pdat3$PC1 + pdat3$PC2 + pdat3$PC3,
             family = binomial())
  pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
  beta <- coef(mod)["gdat2[, j]"]
  ret <- c(j, beta, pval)
  names(ret) <- c("SNP_num", "beta", "pval")
  ret
}
stopCluster(cl)
df.genpc13 <- data.frame("chr"=legend[output.genpc13[,"SNP_num"],"chr"],
                      "bp"=legend[output.genpc13[,"SNP_num"],"position"],
                      "p"=output.genpc13[,"pval"])

## Gender + PC1-5 ##
cl <- makeCluster(4)
registerDoParallel(cl)

output.genpc15 <- foreach(j = 1:nSNPsRemaining, .combine = rbind) %dopar% {
  mod <- glm(pdat3$Y ~ gdat2[, j] + pdat3$gender + pdat3$PC1 + pdat3$PC2 + pdat3$PC3 + pdat3$PC4 + pdat3$PC5,
             family = binomial())
  pval <- summary(mod)$coeff["gdat2[, j]", "Pr(>|z|)"]
  beta <- coef(mod)["gdat2[, j]"]
  ret <- c(j, beta, pval)
  names(ret) <- c("SNP_num", "beta", "pval")
  ret
}
stopCluster(cl)
df.genpc15 <- data.frame("chr"=legend[output.genpc15[,"SNP_num"],"chr"],
                      "bp"=legend[output.genpc15[,"SNP_num"],"position"],
                      "p"=output.genpc15[,"pval"])

# for (j in 1:nSNPsRemaining) {
#     if (j %in% floor(quantile(seq(nSNPsRemaining),seq(0,1,0.1)))) {
#         print(paste("Progress:", names(which(j==floor(quantile(seq(nSNPsRemaining),seq(0,1,0.1)))))), quote=F)
#     }
#     noadj.mod <- summary(glm(pdat$Y~gdat[,j]+pdat$gender,family=binomial()))$coeff
#
# }

#####
#Finding Significant SNPs
#####

alpha = 0.05/nSNPsRemaining
SigPos.noadj = df.noadj$position[df.noadj$p<alpha]
SigPos.gender = df.gender$position[df.gender$p<alpha]

SigPos.anc = df.anc$position[df.anc$p<alpha]
SigSNPs.anc = subset(legend,position %in% SigPos.anc)
PVals.anc = data.frame(pos=df.anc$position[df.anc$p<alpha],
                       pvals=df.anc$p[df.anc$p<alpha])
PVals.anc$pval7 = PVals.anc$pvals * 10^(7)
PVals.anc$beta = output.anc[output.anc[,"pval"]<alpha,"beta"]

SigPos.genanc = df.genanc$position[df.genanc$p<alpha]
SigSNPs.genanc = subset(legend,position %in% SigPos.genanc)
PVals.genanc = data.frame(pos=df.genanc$position[df.genanc$p<alpha],
                       pvals=df.genanc$p[df.genanc$p<alpha])
PVals.genanc$pval7 = PVals.genanc$pvals * 10^(7)
PVals.genanc$beta = output.genanc[output.genanc[,"pval"]<alpha,"beta"]

SigPos.pc13 = df.pc13$position[df.pc13$p<alpha]
SigSNPs.pc13 = subset(legend,position %in% SigPos.pc13)
SigPos.pc15 = df.pc15$position[df.pc15$p<alpha]
SigSNPs.pc15 = subset(legend,position %in% SigPos.pc15)
SigPos.genpc13 = df.genpc13$position[df.genpc13$p<alpha]
SigSNPs.genpc13 = subset(legend,position %in% SigPos.genpc13)
SigPos.genpc15 = df.genpc15$position[df.genpc15$p<alpha]
SigSNPs.genpc15 = subset(legend,position %in% SigPos.genpc15)

SigPos.anc.fac = df.anc.fac$position[df.anc.fac$p<alpha]
SigSNPs.anc.fac = subset(legend,position %in% SigPos.anc.fac)

#######
#Plots#
#######


#### Manhattan Plots ####

## Main Model ##
manhattan(df.anc,chr="chr",bp="position",p="p",
          main="Adjusted for Ancestry",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-Anc.png')
dev.off()

## Under Adjusted Models ##
manhattan(df.noadj,chr="chr",bp="position",p="p",
          main="Unadjusted",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-NoAdj.png')
dev.off()
manhattan(df.gender,chr="chr",bp="position",p="p",
          main="Adjusted for Gender",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-Gen.png')
dev.off()

## Over Adjusted Models ##
manhattan(df.genanc,chr="chr",bp="position",p="p",
          main="Adjusted for Gender and Ancestry",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-GenAnc.png')
dev.off()
manhattan(df.pc13,chr="chr",bp="position",p="p",
          main="Adjusted for PCs 1-3",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-PC13.png')
dev.off()
manhattan(df.pc15,chr="chr",bp="position",p="p",
          main="Adjusted for PCs 1-5",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-PC15.png')
dev.off()
manhattan(df.genpc13,chr="chr",bp="position",p="p",
          main="Adjusted for Gender and PCs 1-3",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-GenPC13.png')
dev.off()
manhattan(df.genpc15,chr="chr",bp="position",p="p",
          main="Adjusted for Gender and PCs 1-5",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-GenPC15.png')
dev.off()
manhattan(df.anc.fac,chr="chr",bp="position",p="p",
          main="Adjusted for Ancestry, No Assumed Inheritance Mode",
          suggestiveline = -log10(.05/nSNPsRemaining))
dev.copy(png,'../Figures/Manh-Anc-Fac.png')
dev.off()


#### Q-Q Plots ####

## Main Model ##
qq(df.anc$p, main="Adjusted for Ancestry")
dev.copy(png,'../Figures/QQ-Anc.png')
dev.off()

## Under Adjusted Models ##
op = par(no.readonly = TRUE)
par(mfrow=c(1,2))
qq(df.noadj$p, main="Unadjusted")
qq(df.gender$p, main="Adjusted for Gender")
par(op)
dev.copy(png,'../Figures/QQ-Under.png')
dev.off()

## Over Adjusted Models ##
qq(df.genanc$p, main="Adj. for Gender, Ancestry")
dev.copy(png,'../Figures/QQ-GenAnc.png')
dev.off()

op = par(no.readonly = TRUE)
par(mfrow=c(2,2))
qq(df.pc13$p, main="Adj. for PCs 1-3")
qq(df.pc15$p, main="Adj. for PCs 1-5")
qq(df.genpc13$p, main="Adj. for Gender and PCs 1-3")
qq(df.genpc15$p, main="Adj. for Gender and PCs 1-5")
par(op)
dev.copy(png,'../Figures/QQ-PCs.png')
dev.off()

qq(df.anc.fac$p, main="Adj. for Ancestry, No Assumed Inheritance Mode")
dev.copy(png,'../Figures/QQ-Anc-Fac.png')
dev.off()

