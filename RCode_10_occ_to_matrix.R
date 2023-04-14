# -------------------------------------------------------------------------- #
#
# This code converts the occurrence sheet to a matrix for the incidence data
#
# Version 13/04/2023
#
# -------------------------------------------------------------------------- #
#
# First check your working directory and set your working directory as appropriate
#
# -------------------------------------------------------------------------- #
#
# Open required libraries for the script
library("reshape")
library("dplyr")

# Read the csv file as my_book
my_book <- read.csv("Data_R/BFO_Meishan.csv")

# Subselect occurrences where the bed assignment is known
# Subselect occurrences where the fossil was identified to species-level
# remove occurrence where the exact bed/subbed is not known
my_book <- my_book[my_book$stg_ass == "TRUE",] 
my_book <- my_book[my_book$sp_assignment == "TRUE",]
my_book <- my_book[my_book$range == "FALSE",]

# Convert occurrence list (my_book) to a matrix.
GCGM.Taxa <- cast(my_book, stg~confirmed_name, fun.aggregate = sum, value = "Count")
row.names(GCGM.Taxa) <- GCGM.Taxa$stg

# There is no real abundance data, so change the counts to 1s (presence)
GCGM.Taxa[GCGM.Taxa > 0.5] <- 1

# Create properties data.frame
GCGM.Prop <- read.csv("stages_Meishan.csv")
row.names(GCGM.Prop) <- GCGM.Prop$stg

# merge properties and taxa datasets
# remove surplus data
GCGM.Genus.Matrix <- merge(GCGM.Taxa, GCGM.Prop, by = 0)
GCGM.Genus.Matrix$stg.x <- NULL
GCGM.Genus.Matrix$stg.y <- NULL

# save incdience matrix as csv file
write.csv(GCGM.Genus.Matrix, "species_matrix_meishan.csv")
