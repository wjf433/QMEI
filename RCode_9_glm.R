# -------------------------------------------------------------------------- #
#
# This code runs the generalized linear models
#
# Version 13/04/2023
#
# -------------------------------------------------------------------------- #
#
# First check your working directory and set your working directory as appropriate
#
# -------------------------------------------------------------------------- #
#Install dependencies, if you do not have this use the install.packages("") command.
library(car)
library(jtools)

# -------------------------------------------------------------------------- #
# Set up the data
#read in the data, proxi=geochemical proxy data
proxi <- read.csv("Data_Imputation/proxies_regression.csv")

#subselect beds 22 to 32
short <- proxi[proxi$Sbottom >= -443,]
short <- short[short$Stop <= 22,]

#subset the geochemical proxies that showed a significant relationship to the diversity measure
#remove na values
data_div <- subset(short, select = c("divRT", "Xccarb", "XO", "XN", "XCd", "XCa"))
data_div <- na.omit(data_div)

# -------------------------------------------------------------------------- #
# Species Richness model

#Run the glm() model, but then check if the model is negatively affected by those variables that we know strongly correlate with each other
#show the model output
#use the car package to investigate the value inflated factors for each parameter
#convert SE to confidence intervals
model <- glm(divRT ~Xccarb+XO+XCa+XCd+XN, data = data_div, family="poisson")
summary(model)
car::vif(model)

#Model 1 XCd dropped because it is highly correlated with XO. XN dropped due to high VIF
model <- glm(divRT ~Xccarb+XO+XCa, data = data_div, family="poisson")
summary(model)
car::vif(model)
summ(model, confint = TRUE, ci.width = .95)

#Model 2 XO dropped because it is highly correlated with XCd.
model <- glm(divRT ~Xccarb+XCd+XCa, data = data_div, family="poisson")
summary(model)
car::vif(model)
summ(model, confint = TRUE, ci.width = .95)


