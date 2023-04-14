#load dependent packages
library("vegan")
library("BiodiversityR")
library("ggsci")
library("ggrepel")
library("ggforce")
library("ggplot2")
library("ggrepel")

#+++Read_In_Data+++#

#Read in species composition data, then subselect beds 22-29a and then subselect the species from the matrix
spe <- read.csv("species_matrix_meishan.csv", row.names=1)
spe <- spe[25:40,]
spe <- spe[,2:603]

#Read in the environmental proxy data, then subselect beds 22-29a and then subselect only proxy data
env <- read.csv("proxies_short_regression.csv", row.names=1, sep=",")
env <- env[25:40,]
env <- env[,1:13]

#+++Partial_dbRDA+++#

#explore the different distance indices
rankindex(env, spe, indices = c("euc", "man", "bray", "jaccard", "kulczynski"),stepacross= FALSE, method =
            "spearman") 

#undertake the initial capscale(partial dbRDA) analysis. 
#can also change the analysis to cca, rda, dbRDA
#can also change the distance measure
dbRDA <- capscale(spe ~ Xccarb+XCorg+XO+XSr+XN+XOs+XFeHr+XFePy+XCd+XS+XLi+XCa+XTOC, env,
               dist="jaccard")

#attach the environmental data, this is required for the contour plot later on.
attach(env)

#ANOVA like permutation test  to assess the significance of constraints.
anova.cca(dbRDA)

#Investigate if multicorrelations bewteen proxies impacted the model performance
vif.cca(dbRDA)

#re-run model dropping XN and XCd because the have high VIFs and strongly correlate with oxygen isotopes 
#and re-investigate the model
dbRDA=capscale(spe ~ Xccarb+XO+XSr+XOs+XFeHr+XFePy+XS+XLi+XCa+XTOC, env,
               dist="jaccard")
vif.cca(dbRDA)
anova.cca(dbRDA)
anova.cca(dbRDA, by="terms", permu=999)

#Now re-run model dropping insignificant factors
#reinvestigate the model
dbRDA=capscale(spe ~ Xccarb+XO+XSr, env,
               dist="jaccard")
vif.cca(dbRDA)
anova.cca(dbRDA)
anova.cca(dbRDA, by="terms", permu=999)

#do the dbRDA plot
plot2 <- ordiplot(dbRDA, choices=c(1,2))

#extract the  sample scores for plotting
sites.long3 <- sites.long(plot2, env.data=env)

#rename the colnames for env, because we will now want to plot the strings
colnames(sites.long3)[1:13]  <- c("δ13C", "δ13Corg","δ18O","87Sr/86Sr", "δ15N", "δ187Os","FeHR/FeT", "FePy/FeHr","δ114/110Cd","δ34S","δ7Li ","δ44/40Ca", "TOC")

#extract the species scores for plotting
spec.envfit <- envfit(plot2, env=spe, permutations=99)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot2, spec.data=spec.data.envfit)
species.long3 <- species.long2[species.long2$r >= 0.6, ]
species.long3$labels <- make.cepnames(species.long3$labels)

#extract the geochemical proxy scores for plotting
vectors.envfit <- envfit(plot2, env=env)
vectors.long3 <- vectorfit.long(vectors.envfit)

#rename the proxys for the plot
vectors.long3$vector <- c("δ13C", "δ13Corg","δ18O","87Sr/86Sr", "δ15N", "δ187Os","FeHR/FeT", "FePy/FeHr","δ114/110Cd","δ34S","δ7Li ","δ44/40Ca", "TOC")

#extract scores for oxygen isotope values
temp.surface <- ordisurf(plot2, y=XO)
temp.grid <- ordisurfgrid.long(temp.surface)

#+++CAPSCALE_Plot+++#

#plot the capscale results
ggplot5 <- ggplot() + 
  geom_contour_filled(data=temp.grid, 
                      aes(x=x, y=y, z=z)) +
  scale_fill_brewer(palette = "Reds", type='div', direction = -1) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(element_text("CAP 1")) +
  ylab(element_text("CAP 2")) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=sites.long3, 
             aes(x=axis1, y=axis2), 
             colour="black", alpha=0.7, size=3) +
  geom_label_repel(data=subset(sites.long3), 
                  aes(x=axis1*1, y=axis2*1, label=labels),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = "white") +
  geom_segment(data=subset(vectors.long3, vector %in% c("δ18O", "δ13C", "87Sr/86Sr")),
               aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
               colour="black", size=2, arrow=arrow()) +
  geom_label_repel(data=subset(vectors.long3, vector %in% c("δ18O", "δ13C", "87Sr/86Sr")), 
                   aes(x=axis1*1, y=axis2*1, label=vector),
                   box.padding   = 0.7, 
                   point.padding = 0.7,
                   color = "black") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#+++Save_Plot+++#

#save the plot
ggsave("dbRDA.pdf", plot = ggplot5, device = "pdf",
       path = NULL, width = 18, height = 15, units = c("cm"),
       dpi = 1000)
