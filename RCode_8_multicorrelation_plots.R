# -------------------------------------------------------------------------- #
#
# The code is to create the multicorrelation plots for each geochemical proxy
#
# -------------------------------------------------------------------------- #
#
# Version 12/04/2023
#
# -------------------------------------------------------------------------- #
#
# First check your working directory and set your working directory as appropriate
#
# -------------------------------------------------------------------------- #
# This code requires the functions within the corrplot package

library("corrplot")

# First check correlation values before using the imputed data
# Read in the average values for each bed
avg <- read.csv("Data_Imputation/proxies_averages.csv")

# Subselect the geochemical proxies
avg <- subset(avg, select = c("Xccarb","XCorg","XO","XSr","XN","XOs","XFeHr", "XFePy", "XCd","XS33","XS34","XLi","XCa", "XTOC", "XZn", "XThU", "OCe"))

# Compute correlation values for all non-missing pairs. Insignificant p-values will be crossed-out
p9 <- cor(avg, method = "pearson", use = "pairwise.complete.obs")
testRes = cor.mtest(p9, conf.level = 0.95)
corrplot.mixed(p9, p.mat = testRes$p, sig.level = 0.05)

# -------------------------------------------------------------------------- #
# First check correlation values using the imputed data, we're using the segmented regression values here

#read in the average values for each bed
reg <- read.csv("Data_Imputation/proxies_regression.csv")

# Subselect the geochemical proxies
reg <- subset(reg, select = c("Xccarb","XCorg","XO","XSr","XN","XOs","XFeHr", "XFePy", "XCd","XS33","XS34","XLi","XCa", "XTOC", "XZn", "XThU", "OCe"))

# Compute correlation values for all non-missing pairs. Insignificant p-values will be crossed-out
p10 <- cor(reg, method = "pearson", use = "pairwise.complete.obs")
testRes = cor.mtest(p10, conf.level = 0.95)
corrplot.mixed(p10, p.mat = testRes$p, sig.level = 0.05)

# -------------------------------------------------------------------------- #
# Repeat, but instead just subselect the beds that will be used in the GLM

reg <- read.csv("Data_Imputation/proxies_regression.csv")

# Subselect the geochemical proxies
reg <- subset(reg, select = c("Xccarb","XCorg","XO","XSr","XN","XOs","XFeHr", "XFePy", "XCd","XS33","XS34","XLi","XCa", "XTOC", "XZn", "XThU", "OCe"))

# Subselect beds 22 - 29a
reg <- reg[25:40,]

# Compute correlation values for all non-missing pairs. Insignificant p-values will be crossed-out
p11 <- cor(reg, method = "pearson", use = "pairwise.complete.obs")
testRes = cor.mtest(p11, conf.level = 0.95)
corrplot.mixed(p11, p.mat = testRes$p, sig.level = 0.05)
