# -------------------------------------------------------------------------- #
#
# This code runs the extinction pulses analyses for each phylum
# Please note that this analysis was run using 90 cores and 600Gb memory,
# The analysis will be slower on personal computer with fewer cores,
# If you computer memory is too small, the code will crash as it cannot store all of the partitions.
#
# Version 22/03/2023
#
# -------------------------------------------------------------------------- #
#
# This code requires the functions to be defined first, which are available from
# RCode_1_Ext_Pulses_Functions.R
#
# -------------------------------------------------------------------------- #
#
# First check your working directory and set your working directory as appropriate
#
# -------------------------------------------------------------------------- #
# Foraminifera

# read in the data
data <- matrix(scan("Data_R/data_foraminifera.txt"), 39, 125, byrow=T)

for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 5            # Set the number of extinction pulses
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)
plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("Figures/foraminifera.pdf",         # File name
    width = 7, height = 6,      # Width and height in inches
    bg = "white",               # Background color
    colormodel = "cmyk",        # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# -------------------------------------------------------------------------- #
# Ostracods

# read in the data
data <- matrix(scan("Data_R/data_arthropoda.txt"), 36, 76, byrow=T)

for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 6            # Set the number of extinction pulses
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)
plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("Figures/arthropoda.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# -------------------------------------------------------------------------- #
# Conodonts

# read in the data
data <- matrix(scan("Data_R/data_chordata.txt"), 59, 49, byrow=T)

for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 5            # Set the number of extinction pulses
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)
plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("Figures/conodonta.pdf",         # File name
    width = 7, height = 6,   # Width and height in inches
    bg = "white",            # Background color
    colormodel = "cmyk",     # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# -------------------------------------------------------------------------- #
# Molluscs

# read in the data
data <- matrix(scan("Data_R/data_mollusca.txt"), 40, 39, byrow=T)

for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 6            # Set the number of extinction pulses
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)
plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("Figures/mollusca.pdf",         # File name
    width = 7, height = 6,  # Width and height in inches
    bg = "white",           # Background color
    colormodel = "cmyk",    # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# -------------------------------------------------------------------------- #
# Brachiopoda

# read in the data
data <- matrix(scan("Data_R/data_brachiopoda.txt"), 34, 33, byrow=T)

for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 6            # Set the number of extinction pulses
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)
plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("Figures/brachiopoda.pdf",         # File name
    width = 7, height = 6,     # Width and height in inches
    bg = "white",              # Background color
    colormodel = "cmyk",       # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()