## Author: Pedro de los Reyes Rodriguez
## Email: pedro.reyes@ibvf.csic.es
## Date: February 2017


##Arguments:
arguments <- commandArgs(trailingOnly=TRUE)

qpcr.table <- arguments[1] ##Table with the qpcr data. COlumns separated by tab. File.txt
main.title <- arguments[2] ##Title of the work (For example: site CORE1 from promoter GBSS)
ylim.bottom <- as.numeric(arguments[3]) ## Bottom Y-axis limit in the plot
ylim.up <- as.numeric(arguments[4]) ## Y-axis upper limit in the plot

# qpcr.data <- read.table(file = "core_pFT.txt", header = TRUE, sep = "\t")
# main.title <- "prueba"
# ylim.bottom <- 0
# ylim.up <- 5


#Extracting the data for INPUT and CHIP samples for each genotype (wild-type and tagged)


qpcr.data <- read.table(file = qpcr.table, header = TRUE, sep = "\t")

col.input <- qpcr.data$Ct[qpcr.data$Genotype == "WT" & qpcr.data$Treatment == "INPUT"]
col.chip <- qpcr.data$Ct[qpcr.data$Genotype == "WT" & qpcr.data$Treatment == "CHIP"]

ctapi.input <- qpcr.data$Ct[qpcr.data$Genotype == "CTAPI" & qpcr.data$Treatment == "INPUT"]
ctapi.chip <- qpcr.data$Ct[qpcr.data$Genotype == "CTAPI" & qpcr.data$Treatment == "CHIP"]




#Applying the ......

i <- 1
enrichment <- c()  

chip <- function(input, chip)
{
  for (i in 1:length(input))
  {
    enrichment[i] <- 10*(2^(input[i]-chip[i]))
    
  }
  
  return(enrichment)
}


ctapi <- chip(ctapi.input, ctapi.chip)
ctapi <- ctapi[!is.na(ctapi)]

col <-chip(col.input,col.chip)
col <- col[!is.na(col)]

ctapi.mean <- mean(ctapi, na.rm = TRUE)
col.mean <- mean(col, na.rm = TRUE)

if (length(ctapi) > 1)
{
  ctapi.sd <- sd(ctapi)
} else if (length(ctapi == 1))
{
  ctapi.sd <- 0
}


if (length(col) > 1)
{
  col.sd <- sd(col)
} else if (length(col == 1))
{
  col.sd <- 0
}


exprs <- c(ctapi.mean, col.mean)
sds <- c(ctapi.sd, col.sd)

# tiff(filename = main.title,width = 10, height = 10, units = "in", res = 600)

#Abro un canal para generar el archivo png
png(file=main.title,
    width     = 10,
    height    = 10,
    units     = "in",
    res       = 600
)

## Genero la imagen 

xpos <- barplot(exprs, names.arg = c("35S:CO CTAPi", "Wild-type"),
                col = c("red", "green"), ylim = c(ylim.bottom,ylim.up), ylab="Enrichment (A.U)",
                main = main.title)


arrows(xpos, exprs-sds, xpos, exprs+sds, angle=90, code=3, length=0.1, lwd = 2)
##Cierro el canal
dev.off()






