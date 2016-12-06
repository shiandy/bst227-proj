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

# hardy-weinberg by ethnicity
hwe_african <- apply(gdat2[pdat2$ancestry == "African", ], 2, hwe_test)
hwe_euro <- apply(gdat2[pdat2$ancestry == "European", ], 2, hwe_test)

hist(hwe_african)
hist(hwe_euro)

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
    geom_point()

ggplot(data = g_pca_df, aes(x = PC2, y = PC3, color = pdat2$ancestry)) +
    geom_point()

plot(svd_res$d, type = "b")



#####################
#Logistic Regression#
#####################


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

# for (j in 1:nSNPsRemaining) {
#     if (j %in% floor(quantile(seq(nSNPsRemaining),seq(0,1,0.1)))) {
#         print(paste("Progress:", names(which(j==floor(quantile(seq(nSNPsRemaining),seq(0,1,0.1)))))), quote=F)
#     }
#     noadj.mod <- summary(glm(pdat$Y~gdat[,j]+pdat$gender,family=binomial()))$coeff
#
# }


#######
#Plots#
#######

manhattan(df.gender,chr="chr",bp="position",p="p", main="Adjusted Covariates: Gender", suggestiveline = -log10(.05/nSNPsRemaining))
qq(df.gender$p, main="Adjusted Covariates: Gender")
