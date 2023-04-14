# -------------------------------------------------------------------------- #
#
# The code is to create the range charts for each species used in Fig. 1
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
# This code requires the functions within the divDyn package
library("divDyn")

# Read in the occurrence and stratigraphic data and subselect the occurrences that will be used in the analysis
fossils <- read.csv("Data_R/BFO_Meishan.csv")
stages <- read.csv("stages_Meishan.csv")


fossils <- fossils[fossils$stg_ass == "TRUE",]       # select only occurrences where the bed is known
fossils <- fossils[fossils$sp_assignment == "TRUE",] # select only the occurrences that are identified to species-level
fossils <- fossils[fossils$range == "FALSE",]        # select only occurrences that are known from a single bed rather than a range of beds

nrow(fossils)                                        # show number of occurrences included in the analysis

# Add the used stratigraphic information to the fossil occurrences, counter allows you to know how well it is doing
for (i in 1: nrow(fossils)){
  x <- fossils[i,"stg"]
  print(i)
  for (j in 1: nrow(stages)){
    y <- stages[j, "stg"]
    if(x==y){
      fossils[i,"Stop"] <- stages[j, "Stop"]
      break
      }
  }
}

# -------------------------------------------------------------------------- #
# Foraminifera

forams <- fossils[fossils$Phylum == "Foraminifera",]  #select only the foraminifera

# create a dataframe that will be used to make the range chart
fl <- fadlad(forams, bin="Stop", tax="confirmed_name")
fl <- fl[fl$duration != 0, ]                           #remove singletons

# create a pdf file to save the plot
pdf("Figures/ranges_foraminifera.pdf",        # File name
    width = 12, height = 8.0354,      # Width and height in inches
    bg = "white",                     # Background color
    colormodel = "cmyk")              # Color model

# create the plot field
plot<- tsplot(stages, shading="Stop", xlim=c(-4390, 773), xlab="cm below Permian/Triassic boundary")
abline(v=22)  # add line at top of bed 29a

ranges(fl, tax="confirmed_name", bin="Stop", labs=F,
       labels.args=list(cex=0.9), occs=F) #plot the ranges of each species
dev.off()  # Close the pdf file

# -------------------------------------------------------------------------- #
# Brachiopoda

brachs <- fossils[fossils$Phylum == "Brachiopoda",] #select only the brachiopods

# create a dataframe that will be used to make the range chart
fl <- fadlad(brachs, bin="Stop", tax="confirmed_name")
fl <- fl[fl$duration != 0, ]                         #remove singletons

# create a pdf file to save the plot
pdf("Figures/ranges_brachiopoda.pdf",         # File name
    width = 12, height = 3.6324,      # Width and height in inches
    bg = "white",                     # Background color
    colormodel = "cmyk")              # Color model

# create the plot field
plot<- tsplot(stages, shading="Stop", xlim=c(-4390, 773), xlab="cm below Permian/Triassic boundary")
abline(v=22)    # add line at top of bed 29a
abline(v=-3268) # add line at top of bed 9

ranges(fl, tax="confirmed_name", bin="Stop", labs=F,
       labels.args=list(cex=0.9), occs=F) #plot the ranges of each species
dev.off() # Close the pdf file

# -------------------------------------------------------------------------- #
# Mollusca

molluscs <- fossils[fossils$Phylum == "Mollusca",] #select only the molluscs

# create a dataframe that will be used to make the range chart
fl <- fadlad(molluscs, bin="Stop", tax="confirmed_name")
fl <- fl[fl$duration != 0, ] #remove singletons

# create a pdf file to save the plot
pdf("Figures/ranges_mollusca.pdf",            # File name
    width = 12, height = 3.7804,      # Width and height in inches
    bg = "white",                     # Background color
    colormodel = "cmyk")              # Color model

# create the plot field
plot<- tsplot(stages, shading="Stop", xlim=c(-4390, 773), xlab="cm below Permian/Triassic boundary")
abline(v=12) # add line at top of bed 28
abline(v=-3268) # add line at top of bed 9

ranges(fl, tax="confirmed_name", bin="Stop", labs=F,
       labels.args=list(cex=0.9), occs=F) #plot the ranges of each species
dev.off() # Close the pdf file

# -------------------------------------------------------------------------- #
# Conodonta

conos <- fossils[fossils$Phylum == "Chordata",] #select only the conodonts

#create a dataframe that will be used to make the range chart
fl <- fadlad(conos, bin="Stop", tax="confirmed_name")
fl <- fl[fl$duration != 0, ] #remove singletons

#create a pdf file to save the plot
pdf("Figures/ranges_conodonta.pdf",           # File name
    width = 12, height = 3.7804,      # Width and height in inches
    bg = "white",                     # Background color
    colormodel = "cmyk")              # Color model

#create the plot field
plot<- tsplot(stages, shading="Stop", xlim=c(-4390, 773), xlab="cm below Permian/Triassic boundary")
abline(v=22) # add line at top of bed 29a
abline(v=-3268) # add line at top of bed 9

ranges(fl, tax="confirmed_name", bin="Stop", labs=F,
       labels.args=list(cex=0.9), occs=F) #plot the ranges of each species
dev.off() # Close the pdf file

# -------------------------------------------------------------------------- #
# Ostracoda and 1 lonely trilobite species

ostracods <- fossils[fossils$Phylum == "Arthropoda",] #select only the ostracods

#create a dataframe that will be used to make the range chart
fl <- fadlad(ostracods, bin="Stop", tax="confirmed_name")
fl <- fl[fl$duration != 0, ] #remove singletons

#create a pdf file to save the plot
pdf("Figures/ranges_arthropoda.pdf",          # File name
    width = 12, height = 5.7044,      # Width and height in inches
    bg = "white",                     # Background color
    colormodel = "cmyk")              # Color model

#create the plot field
plot<- tsplot(stages, shading="Stop", xlim=c(-4390, 773), xlab="cm below Permian/Triassic boundary")
abline(v=-29) # add line at top of bed 24d
abline(v=-171) # add line at top of bed 23a
abline(v=-2595) # add line at top of bed 12

ranges(fl, tax="confirmed_name", bin="Stop", labs=F,
       labels.args=list(cex=0.9), occs=F) #plot the ranges of each  species
dev.off() # Close the pdf file

# -------------------------------------------------------------------------- #
# Everything else

# select only the remaining groups
other <- fossils[fossils$Phylum == "Problematicum" | fossils$Phylum == "Chlorophyta" | fossils$Phylum == "Rhodophyta"| fossils$Phylum == "Bryozoa" | fossils$Phylum == "Cnidaria" | fossils$Phylum == "Radiolaria",]

# create a dataframe that will be used to make the range chart
fl <- fadlad(other, bin="Stop", tax="confirmed_name")
fl <- fl[fl$duration != 0, ] #remove singletons

# create a pdf file to save the plot
pdf("Figures/ranges_other.pdf",               # File name
    width = 12, height = 2.6334,      # Width and height in inches
    bg = "white",                     # Background color
    colormodel = "cmyk")              # Color model

# create the plot field
plot<- tsplot(stages, shading="Stop", xlim=c(-4390, 773), xlab="cm below Permian/Triassic boundary")

ranges(fl, tax="confirmed_name", bin="Stop", labs=F,
       labels.args=list(cex=0.9), occs=F) #plot the ranges of each  species
dev.off() #Close the pdf file
