## Author: Pedro de los Reyes Rodriguez
## Email: pedro.reyes@ibvf.csic.es
## Date: February 2017

##This script works with arguments in order to submit it using Rscript command.
##If you want to run it in Rstudio, don't use arguments. 

#Arguments:
arguments <- commandArgs(trailingOnly=TRUE)

qpcr.table <- arguments[1] ##Table with the qpcr data. Columns separated by tab. File.txt
main.title <- arguments[2] ##Title of the work (For example: site Gbox from GI promoter)
ylim.bottom <- as.numeric(arguments[3]) ## Bottom Y-axis limit in the plot
ylim.up <- as.numeric(arguments[4]) ## Y-axis upper limit in the plot

# qpcr.data <- read.table(file = "pGI_gbox.txt", header = TRUE, sep = "\t")
# main.title <- "GI promoter"
# ylim.bottom <- 0
# ylim.up <- 0.05


## Extracting the data for INPUT and CHIP samples for each genotype (wild-type and tagged)

qpcr.data <- read.table(file = qpcr.data, header = TRUE, sep = "\t")

gfphy5.input <- qpcr.data$Ct[qpcr.data$Genotype == "gfp_hy5" & qpcr.data$Treatment == "INPUT"]
gfphy5.chip <- qpcr.data$Ct[qpcr.data$Genotype == "gfp_hy5" & qpcr.data$Treatment == "CHIP"]

col.input <- qpcr.data$Ct[qpcr.data$Genotype == "col" & qpcr.data$Treatment == "INPUT"]
col.chip <- qpcr.data$Ct[qpcr.data$Genotype == "col" & qpcr.data$Treatment == "CHIP"]

## Adjusting the input to 100% ##
## In this case, the starting input fraction is 10%, then a dilution factor of 100 
## or 3.321928 cycles (log2 of 100), is substracted from the Ct value of diluted input.

gfphy5.input <- gfphy5.input - log2(10)
col.input <- col.input - log2(10)


## Applying the percent of input method for each technical replicate in QPCR

i <- 1
percents <- c()  

percent.input <- function(input, chip)
{
  for (i in 1:length(input))
  {
    percents[i] <- 100*(2^(input[i]-chip[i])) ##I multiply by 100 just to calculate the percentage
  }
  
  return(percents)
}


gfphy5 <- percent.input(gfphy5.input, gfphy5.chip)
gfphy5 <- gfphy5[!is.na(gfphy5)]

col <-percent.input(col.input,col.chip)
col <- col[!is.na(col)]

## Calculate the mean and the sd
gfphy5.mean <- mean(gfphy5, na.rm = TRUE)
col.mean <- mean(col, na.rm = TRUE)

if (length(gfphy5) > 1)
{
  gfphy5.sd <- sd(gfphy5)
} else if (length(gfphy5 == 1))
{
  gfphy5.sd <- 0
}


if (length(col) > 1)
{
  col.sd <- sd(col)
} else if (length(col == 1))
{
  col.sd <- 0
}


exprs <- c(gfphy5.mean, col.mean)
sds <- c(gfphy5.sd, col.sd)

## Barplot

# tiff(filename = main.title,width = 10, height = 10, units = "in", res = 600)

#Generating png file
png(file=main.title,
    width     = 10,
    height    = 10,
    units     = "in",
    res       = 600
)

## Plot

xpos <- barplot(exprs, names.arg = c("gfphy5", "col"),
                col = c("red", "green"), ylim = c(ylim.bottom,ylim.up), ylab="% of Input",
                main = main.title)


arrows(xpos, exprs-sds, xpos, exprs+sds, angle=90, code=3, length=0.1, lwd = 2)
##Close
dev.off()






