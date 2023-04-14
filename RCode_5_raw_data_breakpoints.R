# -------------------------------------------------------------------------- #
#
# The code is to create the breakpoints for species richnes and geochemical proxies
#
# -------------------------------------------------------------------------- #
#
# Version 22/03/2023
#
# -------------------------------------------------------------------------- #
#
# First check your working directory and set your working directory as appropriate
#
# -------------------------------------------------------------------------- #
# This code requires the functions within the ggplot2 and segmented packages
library(ggplot2)
library("segmented")
# -------------------------------------------------------------------------- #
# Species richness changes

# read in the data
env <- read.csv("Data_R/divDyn_meishan.csv")

#subselect the useful columns
env_1 <- subset(env, select = c("Stop","divRT")) #subselect the useful columns
env_1 <- env_1[env_1$Stop >-4389.9,] #subselect the useful part of the section
env_1 <- env_1[env_1$Stop <773.1,]

#determine the number of breakpoints following https://rdrr.io/cran/segmented/man/selgmented.html
out.lm<-lm(divRT~Stop,data=env_1)
os <-selgmented(out.lm, Kmax=10, type="bic") #BIC-based selection

#Do the breakpoint analysis
fit <- lm(divRT ~ Stop, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~Stop, psi=c(-3635, -629, -236, -16, 24))
summary(segmented.fit)
fit <- numeric(length(env_1$Stop)) * NA
fit[complete.cases(rowSums(cbind(env_1$divRT, env_1$divRT)))] <- broken.line(segmented.fit)$fit

#plot the results
data1 <- data.frame(Stop = env_1$Stop, divRT = env_1$divRT, fit = fit)
p1<- ggplot(data=env_1, aes(x=Stop, y=divRT, group=1)) + 
  geom_line(colour="#CCCCCC", linetype="dashed", size=.5) + 
  geom_point(colour="#CCCCCC", size=2, shape=21, fill="white") +
  geom_line(aes(x = Stop, y = fit), size = .75, color = "#999999") +
  xlab("") + ylab("Species richness") +
  scale_x_reverse(limits = c(773,  -560)) +
  scale_y_continuous(breaks=c(0,50,100,150, 200, 250))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_diversity_short.pdf", plot = p1, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Oxygen isotopes

# read in the data
env <- read.csv("Data_R/O_18.csv")
env_1 <- subset(env, select = c("n_height","XO")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints following https://rdrr.io/cran/segmented/man/selgmented.html
out.lm<-lm(XO~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=10, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XO ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-3970, -842, -19, -3))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XO, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XO = env_1$XO, fit = fit)
p4<- ggplot(data=env_1, aes(x=n_height, y=XO, group=1)) + 
  geom_line(colour="#FF0000", linetype="dashed", size=.5) + 
  geom_point(colour="#FF0000", size=2, shape=21, fill="#FFFFFF") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#CC0000") +
  xlab("") + ylab("δ18Oapatite (‰, VSMOW)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_O_isotope.pdf", plot = p4, device = "pdf",
  path = NULL, width = 12.6, height = 4, units = c("cm"),
  dpi = 1000)

# -------------------------------------------------------------------------- #
# Carbonate carbon isotopes

# read in the data
env <- read.csv("Data_R/C_carb_13.csv")
env_1 <- subset(env, select = c("n_height","Xccarb")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(Xccarb~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=10, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(Xccarb ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-3000, -2500, -2000, -1500, -1000, 0))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XCcarb, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XCcarb = env_1$Xccarb, fit = fit)
p1<- ggplot(data=env_1, aes(x=n_height, y=Xccarb, group=1)) + 
  geom_line(colour="#339900", linetype="dashed", size=.5) + 
  geom_point(colour="#339900", size=2, shape=21, fill="#FFFFFF") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#336600") +
  xlab("") + ylab("δ13Ccarb (‰, VPDB)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Ccarb_isotope.pdf", plot = p1, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Organic carbon isotopes

# read in the data
env <- read.csv("Data_R/C_org_13.csv")
env_1 <- subset(env, select = c("n_height","XCorg")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XCorg~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=3, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XCorg ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-500, -4))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XCorg, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
p2<- ggplot(data=env_1, aes(x=n_height, y=XCorg, group=1)) + 
  geom_line(colour="#666666", linetype="dashed", size=.5) +
  geom_point(colour="#666666", size=2, shape=21, fill="#FFFFFF") +
  xlab("") + ylab("δ13Corg (‰, VPDB)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Corg_isotope.pdf", plot = p2, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Strontium isotopes

# read in the data
env <- read.csv("Data_R/Sr_87.csv")
env_1 <- subset(env, select = c("n_height","XSr")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XSr~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XSr ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(172))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XSr, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XSr = env_1$XSr, fit = fit)
p5<- ggplot(data=env_1, aes(x=n_height, y=XSr, group=1)) + 
  geom_line(colour="#0066FF", linetype="dashed", size=.5) + 
  geom_point(colour="#0066FF", size=2, shape=21, fill="#FFFFFF") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#0033FF") +
  xlab("") + ylab("87Sr/86Srapatite (‰)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Sr_isotope.pdf", plot = p5, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Total organic Carbon

# read in the data
env <- read.csv("Data_R/TOC.csv")
env_1 <- subset(env, select = c("n_height","XTOC")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XTOC~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XTOC ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-428))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XSr, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
p3<- ggplot(data=env_1, aes(x=n_height, y=XTOC, group=1)) + 
  geom_line(colour="#000000", linetype="dashed", size=.5) + 
  geom_point(colour="#000000", size=2, shape=21, fill="#FFFFFF") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#000000") +
  xlab("") + ylab("TOC(%)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_TOC.pdf", plot = p3, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Nitrogen isotopes

# read in the data
env <- read.csv("Data_R/N_15.csv")
env_1 <- subset(env, select = c("n_height","XN")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XN~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XN ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-20, -4))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XN, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XSr = env_1$XN, fit = fit)
p6<- ggplot(data=env_1, aes(x=n_height, y=XN, group=1)) + 
  geom_line(colour="#33FF99", linetype="dashed", size=.5) + 
  geom_point(colour="#33FF99", size=2, shape=21, fill="#FFFFFF") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#009966") +
  xlab("") + ylab("δ15N (‰)") +
  scale_x_reverse(limits = c(773, -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_N_isotope.pdf", plot = p6, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Calcium isotopes

# read in the data
env <- read.csv("Data_R/Ca_44.csv")
env_1 <- subset(env, select = c("n_height","XCa")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XCa~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XCa ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(0))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XCa, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XCa = env_1$XCa, fit = fit)
p7<- ggplot(data=env_1, aes(x=n_height, y=XCa, group=1)) + 
  geom_line(colour="#FFCC66", linetype="dashed", size=.5) + 
  geom_point(colour="#FFCC66", size=2, shape=21, fill="white") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#FF9900") +
  xlab("") + ylab("δ44/40Caapatite (‰ BSE)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Ca_isotope.pdf", plot = p7, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Cadmium isotopes

# read in the data
env <- read.csv("Data_R/Cd_114.csv")
env_1 <- subset(env, select = c("n_height","XCd")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XCd~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XCd ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-23, 4))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XCd, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XCd = env_1$XCd, fit = fit)
p8<- ggplot(data=env_1, aes(x=n_height, y=XCd, group=1)) + 
  geom_line(colour="#CCCCCC", linetype="dashed", size=.5) + 
  geom_point(colour="#CCCCCC", size=2, shape=21, fill="white") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#999999") +
  xlab("") + ylab("δ114/110Cd (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Cd_isotope.pdf", plot = p8, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Lithium isotopes

# read in the data
env <- read.csv("Data_R/Li_7.csv")
env_1 <- subset(env, select = c("n_height","XLi")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XLi~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XLi ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(172))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XLi, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XLi = env_1$XLi, fit = fit)
p9<- ggplot(data=env_1, aes(x=n_height, y=XLi, group=1)) + 
  geom_line(colour="#FF99FF", linetype="dashed", size=.5) + 
  geom_point(colour="#FF99FF", size=2, shape=21, fill="white") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#FF66FF") +
  xlab("") + ylab("δ7Li (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Li_isotope.pdf", plot = p9, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Sulphur isotopes

# read in the data
env <- read.csv("Data_R/S_33_34.csv")
env_1 <- subset(env, select = c("n_height","XS34")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XS34~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# plot the results
p10<- ggplot(data=env_1, aes(x=n_height, y=XS, group=1)) + 
  geom_line(colour="#FF9900", linetype="dashed", size=.5) + 
  geom_point(colour="#FF9900", size=2, shape=21, fill="white") +
  xlab("") + ylab("δ34S (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_S_34_isotope.pdf", plot = p10, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# read in the data
env_1 <- subset(env, select = c("n_height","XS33")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XS33~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# plot the results
p10<- ggplot(data=env_1, aes(x=n_height, y=XS, group=1)) + 
  geom_line(colour="#FF9900", linetype="dashed", size=.5) + 
  geom_point(colour="#FF9900", size=2, shape=21, fill="white") +
  xlab("") + ylab("δ34S (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_S_33_isotope.pdf", plot = p10, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Osmium isotopes

# read in the data
env <- read.csv("Data_R/Os_187.csv")
env_1 <- subset(env, select = c("n_height","XOs")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints.
out.lm<-lm(XOs~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# plot the results
p11<- ggplot(data=env_1, aes(x=n_height, y=XOs, group=1)) + 
  geom_line(colour="gray", linetype="dashed", size=.5) + 
  geom_point(colour="gray", size=2, shape=21, fill="#CCFFFF") +
  xlab("") + ylab("δ187Os (‰)") +
  scale_x_reverse(limits = c(773, -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Os_isotope.pdf", plot = p11, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Mercury/TOC ratio

# read in the data
env <- read.csv("Data_R/Hg_TOC.csv")
env_1 <- subset(env, select = c("n_height","XHg")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XHg~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# plot the results
p12<- ggplot(data=env_1, aes(x=n_height, y=XHg, group=1)) + 
  geom_line(colour="black", linetype="dashed", size=.5) + 
  geom_point(colour="black", size=2, shape=21, fill="black") +
  xlab("") + ylab("Hg/TOC") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Hg.pdf", plot = p12, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Iron speciation

# read in the data
env <- read.csv("Data_R/Fe_sp.csv")
env_1 <- subset(env, select = c("n_height","XFeHr")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XFeHr~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XFeHr ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-103, 124))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XLi, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
p13<- ggplot(data=env_1, aes(x=n_height, y=XFeHr, group=1)) + 
  geom_line(colour="#993300", linetype="dashed", size=.5) + 
  geom_point(colour="#993300", size=2, shape=21, fill="white") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "brown") +
  xlab("") + ylab("FeHR/FeT") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_FeHr_isotope.eps", plot = p13, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# read in the data
env_1 <- subset(env, select = c("n_height","XFePy")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XFePy~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# plot the results
p14<- ggplot(data=env_1, aes(x=n_height, y=XFePy, group=1)) + 
  geom_line(colour="#993300", linetype="dashed", size=.5) + 
  geom_point(colour="#993300", size=2, shape=21, fill="white") +
  geom_smooth(method=lm, se=FALSE, size = .75, color = "brown") +
  xlab("") + ylab("FePy/FeHr") +
  scale_x_reverse(limits = c(773, -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_FePy_isotope.eps", plot = p14, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Zinc isotopes

# read in the data
env <- read.csv("Data_R/Zn_66.csv")
env_1 <- subset(env, select = c("n_height","XZn")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XZn~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XZn ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-5, 37, 116))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XCa, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XZn = env_1$XZn, fit = fit)
p7<- ggplot(data=env_1, aes(x=n_height, y=XZn, group=1)) + 
  geom_line(colour="#FFCC66", linetype="dashed", size=.5) + 
  geom_point(colour="#FFCC66", size=2, shape=21, fill="white") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#FF9900") +
  xlab("") + ylab("δ66Zncarb (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Zn_isotope.pdf", plot = p7, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Cerium anomalies

# read in the data
env <- read.csv("Data_R/Ce_anom.csv")
env_1 <- subset(env, select = c("n_height","OCe")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(OCe~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# plot the results
p5<- ggplot(data=env_1, aes(x=n_height, y=OCe, group=1)) + 
  geom_line(colour="#0066FF", linetype="dashed", size=.5) + 
  geom_point(colour="#0066FF", size=2, shape=21, fill="#FFFFFF") +
  xlab("") + ylab("Conodont OCeapatite") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_Ce_anom.pdf", plot = p5, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Thorium/Uranium ratios

# read in the data
env <- read.csv("Data_R/Th_U.csv")
env_1 <- subset(env, select = c("n_height","XThU")) #subselect the useful columns
env_1 <- env_1[env_1$n_height >-4390,] #subselect the useful part of the section
env_1 <- env_1[env_1$n_height <773,]

# determine the number of breakpoints
out.lm<-lm(XThU~n_height,data=env_1)
os <-selgmented(out.lm, Kmax=5, type="bic") #BIC-based selection

# do the breakpoint analysis
fit <- lm(XThU ~ n_height, data=env_1)
segmented.fit <- segmented(fit, seg.Z = ~n_height, psi=c(-79))
summary(segmented.fit)
fit <- numeric(length(env_1$n_height)) * NA
fit[complete.cases(rowSums(cbind(env_1$XSr, env_1$n_height)))] <- broken.line(segmented.fit)$fit

# plot the results
data1 <- data.frame(n_height = env_1$n_height, XThU = env_1$XThU, fit = fit)
p5<- ggplot(data=env_1, aes(x=n_height, y=XThU, group=1)) + 
  geom_line(colour="#0066FF", linetype="dashed", size=.5) + 
  geom_point(colour="#0066FF", size=2, shape=21, fill="#FFFFFF") +
  geom_line(aes(x = n_height, y = fit), size = .75, color = "#0033FF") +
  xlab("") + ylab("Th/Uapatite") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# save the results
ggsave("Figures/breakpoint_ThU_isotope.pdf", plot = p5, device = "pdf",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)
