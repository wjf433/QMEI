# -------------------------------------------------------------------------- #
#
# The code is to create the aggregated and imputed values
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
# Attach bed numbers and calculate averages to raw data sheets

#Read in the geochemical and stratigraphic data
bed <- read.csv("stages_Meishan.csv")

# This creates a list of each raw data file
proxies <- c("Data_R/TOC.csv","Data_R/Th_U.csv", "Data_R/Zn_66.csv", "Data_R/S_33_34.csv", "Data_R/S_33_34.csv","Data_R/Ce_anom.csv", "Data_R/Sr_87.csv",
             "Data_R/Os_187.csv", "Data_R/O_18.csv", "Data_R/N_15.csv", "Data_R/Li_7.csv", "Data_R/Hg_TOC.csv", "Data_R/Fe_sp.csv", "Data_R/Fe_sp.csv", 
             "Data_R/Cd_114.csv", "Data_R/Ca_44.csv", "Data_R/C_org_13.csv", "Data_R/C_carb_13.csv")

# This creates the name of the column of interest for each proxy
pnam <- c("XTOC","XThU", "XZn", "XS33", "XS34", "OCe", "XSr",
          "XOs", "XO", "XN", "XLi", "XHg", "XFeHr", "XFePy", 
          "XCd", "XCa", "XCorg", "Xccarb")

# create a dataframe to save a value for each bed in the section
total <- data.frame("Group.1" = c("0","1","1u","2","3","4a","4b","5","6","7","8","9","10","11","12","13a","13b","14","15","16","17","18","19",
                                  "20","21","22","23a","23b","24a","24b","24c","24d","24e","25","26","27a","27b","27c","27d","28","29a","29b","30","31","32", 
                                  "33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50+"))

for (i in 1: length(proxies)){
  
  env <- read.csv(proxies[[i]])
  env <- subset(env, select = c("n_height",pnam[[i]])) # subselect the useful columns
  
  #Now add the bed numbers to the dataframe
  env$Bed <- 0 
  for (i in 1: nrow(env)){
    
    height <- env[i,"n_height"]
    
    for (j in 1: nrow(bed)){
      
      bas <- bed[j, "Sbottom"]
      top <- bed[j, "Stop"]
      
      if (height >= bas && height <= top){
        
        env[i, "Bed"] <- bed[j, "bed"]
        
      }
    }
  }
  
# this averages the values for each bed
env_agg <- aggregate(env, by=list(env$Bed), FUN=mean)
env_agg <- env_agg[order(env_agg$n_height),] #order the dataframe occording to bed height
env_agg$n_height <- NULL
env_agg$Bed <- NULL

total <- merge(total,env_agg,by="Group.1",all=T)
  
}

# this attaches the stratigraphic data and puts the bed numbers in order
names(total)[names(total) == "Group.1"] <- "bed"
total <- merge(total, bed,by="bed",all=T)
total <- total[order(total$Smid),]
rownames(total) <- total$bed #change the rownames

# Now we have a dataframe with the average values for each bed and NAs where there ar eno values for that bed
# save the results
write.csv(total, "Data_Imputation/proxies_averages.csv")

# -------------------------------------------------------------------------- #
#
# Imputation using 2 point interpolation
#
# -------------------------------------------------------------------------- #
# This code requires the functions within the dplyr, forcats and zoo packages
library(dplyr)
library(zoo)
library(forcats)

# prepare a new datframe and the data
data <- total[,1:19]
data <- data %>% 
  arrange(fct_relevel(bed, "0","1","1u","2","3","4a","4b","5","6","7","8","9","10","11","12","13a","13b","14","15","16","17","18","19",
                      "20","21","22","23a","23b","24a","24b","24c","24d","24e","25","26","27a","27b","27c","27d","28","29a","29b","30","31","32", 
                      "33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50+"))
data <- data[1:62,]

#this creates 2-point interpolation for each NA
avg <- data %>%
  mutate(X = na.approx(data))
avg <- avg[,c(1,20)]
avg <- avg[6:47,] #subselects the dataframe to the beds of interest (beds 4a to 34)
  
# save the results
write.csv(avg, "Data_Imputation/proxies_interpolated.csv")

# -------------------------------------------------------------------------- #
#
# Imputation using segmented regression lines
#
# -------------------------------------------------------------------------- #
# Raw data and average points give different values (naturally) so we make the regression lines from the raw data
# this code is requires more manual input and changes, because you need to enter the segmented regression values

#read in the proxy dataset of interest
# file names are: "Data_R/TOC.csv","Data_R/Th_U.csv", "Data_R/Zn_66.csv", "Data_R/S_33_34.csv", "Data_R/S_33_34.csv","Data_R/Ce_anom.csv", "Data_R/Sr_87.csv",
#"Data_R/Os_187.csv", "Data_R/O_18.csv", "Data_R/N_15.csv", "Data_R/Li_7.csv", "Data_R/Hg_TOC.csv", "Data_R/Fe_sp.csv", "Data_R/Fe_sp.csv", 
#"Data_R/Cd_114.csv", "Data_R/Ca_44.csv", "Data_R/C_org_13.csv", "Data_R/C_carb_13.csv"
env <- read.csv("Data_R/Fe_sp.csv")

#subselect the useful columns, in this example highly reactive iron
env_1 <- subset(env, select = c("n_height","XFeHr"))
env_1 <-na.omit(env_1) # remove NAs, if any

# subselect breakpoints
# using the breakpoint anaylsis enter the breakpoints for the regression line you want to get the interpolated value from
env_1 <- env_1[env_1$n_height > -70.5,]
env_1 <- env_1[env_1$n_height > 98,]

#fit linear regression model using data frame
model <- lm(XFeHr ~ n_height, data = env_1)
summary(model)

#interpolate using the regression model
# enter the bed heights of the missing values you want to calculate
var <- data.frame(n_height = c(-90.5, -76.5))
predict(model, newdata = var) # returns values to add to your spreadsheet.

# -------------------------------------------------------------------------- #
#
# Imputation using predictive mean matching
#
# -------------------------------------------------------------------------- #
# This code requires the functions within the mice package
library("mice")

# prepare a new datframe and the data
data <- total[,1:19]
data <- data %>% 
  arrange(fct_relevel(bed, "0","1","1u","2","3","4a","4b","5","6","7","8","9","10","11","12","13a","13b","14","15","16","17","18","19",
                      "20","21","22","23a","23b","24a","24b","24c","24d","24e","25","26","27a","27b","27c","27d","28","29a","29b","30","31","32", 
                      "33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50+"))
data <- data[1:62,]

# now the mice function creates the pmm values and then we aggregate it into a new dataframe
imputed_Data <- mice(data, m=5, maxit = 50, method = 'pmm', seed = 500)
imp_env <- complete(imputed_Data,"long")
imp_agg <- aggregate(imp_env, by=list(imp_env$.id), FUN=mean)
imp_agg$.imp <- NULL
imp_agg$.id <- NULL
imp_agg$Group.1 <- NULL
imp_agg$bed <- data$bed # put th ebed names back

# save the results
write.csv(imp_agg, "Data_Imputation/proxies_mice.csv")
