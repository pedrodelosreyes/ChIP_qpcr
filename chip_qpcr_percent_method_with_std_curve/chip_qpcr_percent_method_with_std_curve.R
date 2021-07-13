## Author: Pedro de los Reyes Rodriguez
## Email: pedro.reyes@ibvf.csic.es
## Date: March 2021

# A powerful way to determine whether your qPCR assay is 
# optimized is to run serial dilutions of a template and 
# use the results to generate a standard curve. 
# The template used for this purpose can be a target with 
# known concentration (e.g., nanograms of genomic DNA 
# or copies of plasmid DNA) or a sample of unknown quantity 
# (e.g., cDNA). The standard curve is constructed by plotting 
# the log of the starting quantity of template 
# (or the dilution factor, for unknownquantities) 
# against the CTvalue obtained during amplification of each dilution. 
# The equation of the linear regression line, 
# along with Pearson’s correlationcoefficient (r) or 
# the coefficient of determination (R2), can then be used 
# to evaluate whether your qPCR assay is optimized.

# Also, using the standard curve for each primers pair, we take
# into account the pcr efficiency for each one. Thus we can 
# compare different amplicons

quantities <- c(0.001,0.0001,0.00001)
dilution.factors <- c(10,100,1000) 
ct.values <- c(18.4140918285678, 21.877853335891, 25.0173875151059)

std.curve <- data.frame(quantities, dilution.factors, ct.values)
head(std.curve)

# build a simple regression model that we can use to predict 
# quantities (x) by establishing a statistically significant 
# linear relationship with Ct values (y)


# Load required packages:
#   
#   tidyverse for data manipulation and visualization
#   ggpubr: creates easily a publication ready-plot

library(tidyverse)
library(ggpubr)


theme_set(theme_pubr())

ggplot(std.curve, aes(y = ct.values, x = log10(quantities))) +
  geom_point() +
  stat_smooth(method = lm)

# compute the correlation coefficient between the two variables
cor(std.curve$ct.values, std.curve$quantities)

model <- lm(ct.values ~ log10(quantities), data = std.curve)
model
#y=-3.305x + 8.550

summary(model)

other.ct.value <- 19.5170410715608
quantitie.of.interest <- (other.ct.value-model$coefficients[1])/model$coefficients[2]
10^quantitie.of.interest 

#hurra



##Real data:

# In this case, my genotypes are co10 (co mutant, control)
# and co overexpressor (35SCO). You can modify the variable
# names with your condition

########-------ACTINA Promoter-------####

qpcr.data <- read.table(file = "act_promoter.csv", sep = ",", header = TRUE)


std1 <- mean(qpcr.data$ct_values[qpcr.data$sample == "std1"])
std10 <- mean(qpcr.data$ct_values[qpcr.data$sample == "std10"])
std100 <- mean(qpcr.data$ct_values[qpcr.data$sample == "std100"])
quantities <- c(0.001,0.0001,0.00001)

std.curve <- data.frame(quantities, c(std1,std10,std100))
colnames(std.curve) <- c("quantities", "ct.values")
head(std.curve)

ggplot(std.curve, aes(y = ct.values, x = log10(quantities))) +
  geom_point() +
  stat_smooth(method = lm)

# compute the correlation coefficient between the two variables
cor(std.curve$ct.values, std.curve$quantities)

#Build the model
model <- lm(ct.values ~ log10(quantities), data = std.curve)
model


#Extract data, predict starting DNA quantity for each technical rep
co10.input <- qpcr.data$ct_values[qpcr.data$sample == "co10_input"]
co10.ip <- qpcr.data$ct_values[qpcr.data$sample == "co10_ip"]
x35sco.input <- qpcr.data$ct_values[qpcr.data$sample == "35sco_input"]
x35sco.ip <- qpcr.data$ct_values[qpcr.data$sample == "35sco_ip"]



co10.input.quant <- 10^((co10.input-model$coefficients[1])/model$coefficients[2])
co10.ip.quant <- 10^((co10.ip-model$coefficients[1])/model$coefficients[2])
x35sco.input.quant <- 10^((x35sco.input-model$coefficients[1])/model$coefficients[2])
x35sco.ip.quant <- 10^((x35sco.ip-model$coefficients[1])/model$coefficients[2])

#Calculate mean and estimate the result/percent/enrichment or whatever
co10.input.mean <- mean(co10.input.quant, na.rm = TRUE)
co10.ip.mean <- mean(co10.ip.quant, na.rm = TRUE)
co10.result <- co10.ip.mean*10/co10.input.mean

x35sco.input.mean <- mean(x35sco.input.quant, na.rm = TRUE)
x35sco.ip.mean <- mean(x35sco.ip.quant, na.rm = TRUE)
x35sco.result <- x35sco.ip.mean*10/x35sco.input.mean


barplot(c(co10.result, x35sco.result))

########-------FT Promoter-------####

qpcr.data <- read.table(file = "FT_promoter_CORE1.csv", sep = ",", header = TRUE)


std1 <- mean(qpcr.data$ct_values[qpcr.data$sample == "std1"])
std10 <- mean(qpcr.data$ct_values[qpcr.data$sample == "std10"])
std100 <- mean(qpcr.data$ct_values[qpcr.data$sample == "std100"])
quantities <- c(0.001,0.0001,0.00001)

std.curve <- data.frame(quantities, c(std1,std10,std100))
colnames(std.curve) <- c("quantities", "ct.values")
head(std.curve)

ggplot(std.curve, aes(y = ct.values, x = log10(quantities))) +
  geom_point() +
  stat_smooth(method = lm)

# compute the correlation coefficient between the two variables
cor(std.curve$ct.values, std.curve$quantities)

#Build the model
model <- lm(ct.values ~ log10(quantities), data = std.curve)
model


#Extract data, predict starting DNA quantity for each technical rep
co10.input <- qpcr.data$ct_values[qpcr.data$sample == "co10_input"]
co10.ip <- qpcr.data$ct_values[qpcr.data$sample == "co10_ip"]
x35sco.input <- qpcr.data$ct_values[qpcr.data$sample == "35sco_input"]
x35sco.ip <- qpcr.data$ct_values[qpcr.data$sample == "35sco_ip"]



co10.input.quant <- 10^((co10.input-model$coefficients[1])/model$coefficients[2])
co10.ip.quant <- 10^((co10.ip-model$coefficients[1])/model$coefficients[2])
x35sco.input.quant <- 10^((x35sco.input-model$coefficients[1])/model$coefficients[2])
x35sco.ip.quant <- 10^((x35sco.ip-model$coefficients[1])/model$coefficients[2])

#Calculate mean and estimate the result/percent/enrichment or whatever
co10.input.mean <- mean(co10.input.quant, na.rm = TRUE)
co10.ip.mean <- mean(co10.ip.quant, na.rm = TRUE)
co10.result <- co10.ip.mean*10/co10.input.mean

x35sco.input.mean <- mean(x35sco.input.quant, na.rm = TRUE)
x35sco.ip.mean <- mean(x35sco.ip.quant, na.rm = TRUE)
x35sco.result <- x35sco.ip.mean*10/x35sco.input.mean


barplot(c(co10.result, x35sco.result))



