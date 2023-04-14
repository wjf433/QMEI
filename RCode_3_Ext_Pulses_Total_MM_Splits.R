# -------------------------------------------------------------------------- #
#
# This code runs the extinction pulses analyses for each split for all taxa that go
# extinct in the Meishan member and Yinkeng Formation.
# 
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
# Split 1

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_1.txt"), 59, 290, byrow=T)

for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("split_1_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_1_MM.csv")

# -------------------------------------------------------------------------- #
# Split 2

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_2.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("split_2_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_2_MM.csv")

# -------------------------------------------------------------------------- #
# Split 3

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_3.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("split_3_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_3_MM.csv")

# -------------------------------------------------------------------------- #
# Split 4

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_4.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure
pdf("split_4_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_4_MM.csv")

# -------------------------------------------------------------------------- #
# Split 5

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_5.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure 
pdf("split_5_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_5_MM.csv")

# -------------------------------------------------------------------------- #
# Split 6

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_6.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure 
pdf("split_6_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_6_MM.csv")

# -------------------------------------------------------------------------- #
# Split 7

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_7.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure 
pdf("split_7_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_7_MM.csv")

# -------------------------------------------------------------------------- #
# Split 8

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_8.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure 
pdf("split_8_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_8_MM.csv")

# -------------------------------------------------------------------------- #
# Split 9

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_9.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure 
pdf("split_9_MM.pdf",         # File name
    width = 7, height = 6,    # Width and height in inches
    bg = "white",             # Background color
    colormodel = "cmyk",      # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_9_MM.csv")

# -------------------------------------------------------------------------- #
# Split 10

# read in the data
data <- matrix(scan("Data_R/data_meishan_all_MM_split_10.txt"), 59, 290, byrow=T)
for (i in 1: ncol(data)){
  data[,i] <- sort(data[,i], na.last=T, decreasing = T)
}

numtaxa <- dim(data)[2]
ord <- order(data[1,])
data <- data[,ord]           # sort data by highest occurrence
par(mfrow=c(1,1))
maxnumpulses <- 4
rangechart(data, NA)

partitionlist <- vector("list", maxnumpulses)
for(pulse in 1:maxnumpulses)  
  partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 

pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)

pulsesCI(pulseout[,dim(pulseout)[2]])

# Save your results as a figure 
pdf("split_10_MM.pdf",         # File name
    width = 7, height = 6,     # Width and height in inches
    bg = "white",              # Background color
    colormodel = "cmyk",       # Color model
    paper = "A4")  

plotpulses(data, pulseout)
pulsesCI(pulseout[,dim(pulseout)[2]])

dev.off()

# Export the results as a .csv file
write.csv(pulseout, "split_10_MM.csv")