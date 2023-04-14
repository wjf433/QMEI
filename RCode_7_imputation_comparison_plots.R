# -------------------------------------------------------------------------- #
#
# The code is to create plots for each proxy to visually compare the values of the different imputed datasets
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
# This code requires the functions within the ggplot2 and dplyr packages
library(ggplot2)
library(dplyr)

# read in the data
com <- read.csv("Data_Imputation/Imputation_comparison.csv")

# -------------------------------------------------------------------------- #
# Carbonate Carbon isotopes
Cc <- com[com$proxy == "XCcarb",] # subselect the proxy to be plotted

# Plot the data
C <- ggplot(data=Cc, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ13Ccarb (‰, VPDB)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

# Save the plot
ggsave("C_imp_comparison.eps", plot = C, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Oxygen isotopes
XO <- com[com$proxy == "XO",]

O <- ggplot(data=XO, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ18Oapatite (‰, VSMOW)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("O_imp_comparison.eps", plot = O, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Strontium isotopes
XSr <- com[com$proxy == "XSr",]

p5<- ggplot(data=XSr, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("87Sr/86Srapatite (‰)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("Sr_imp_comparison.eps", plot = p5, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Nitrogen isotopes
XN <- com[com$proxy == "XN",]         

p6<- ggplot(data=XN, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ15N (‰)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("N_imp_comparison.eps", plot = p6, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Organic carbon isotopes
XCo <- com[com$proxy == "XCorg",]    

p2<- ggplot(data=XCo, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ13Corg (‰, VPDB)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("Co_imp_comparison.eps", plot = p2, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# OTotal organic carbon
XTOC <- com[com$proxy == "XTOC",]

p5<- ggplot(data=XTOC, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("TOC(%)") +
  scale_x_reverse(limits = c(773,  -4390)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("TOC_imp_comparison.eps", plot = p5, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Calcium isotopes
XCa <- com[com$proxy == "XCa",]

s1<- ggplot(data=XCa, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ44/40Caapatite (‰ BSE)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("Ca_imp_comparison.eps", plot = s1, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Osmium isotopes
XOs <- com[com$proxy == "XOs",]

s2<- ggplot(data=XOs, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ187Os (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("Os_imp_comparison.eps", plot = s2, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Iron Speciation
XFeHr <- com[com$proxy == "XFeHr",]

s3<- ggplot(data=XFeHr, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("FeHR/FeT") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("FeHr_imp_comparison.eps", plot = s3, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

XFePy <- com[com$proxy == "XFePy",]

s4<- ggplot(data=XFePy, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("FePy/FeHr") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("FePy_imp_comparison.eps", plot = s4, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Sulphur isotopes

XS <- com[com$proxy == "XS",]

s5<- ggplot(data=XS, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ34S (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("S_imp_comparison.eps", plot = s5, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Cadmium isotopes
XCd <- com[com$proxy == "XCd",]

s6<- ggplot(data=XCd, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ114/110Cd (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("Cd_imp_comparison.eps", plot = s6, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)

# -------------------------------------------------------------------------- #
# Lithium isotopes
XLi <- com[com$proxy == "XLi",]

s7<- ggplot(data=XLi, aes(x=n_height, y=value, col=type)) + 
  geom_point(aes(shape=type, color=type), size=2)+
  scale_shape_manual(values=c(15, 17, 21, 18))+
  scale_color_manual(values=c('#66CC33','#E69F00', '#56B4E9', "#9900CC"))+
  xlab("") + ylab("δ7Li (‰)") +
  scale_x_reverse(limits = c(773,  -560)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90))

ggsave("Li_imp_comparison.eps", plot = s7, device = "eps",
       path = NULL, width = 12.6, height = 4, units = c("cm"),
       dpi = 1000)
