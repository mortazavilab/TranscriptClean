# Dana Wyman
# 12/29/2017
# This script generates a PDF report based on a TranscriptClean run
# Designed for use with R version R/3.3.2

# Read input arguments
args = commandArgs(trailingOnly = TRUE)
logFile = args[1]
prefix = args[2]
reportFile = paste(prefix, "report.pdf", sep="_")


# Check for packages and install if not found
#print("Checking for R packages...")
#all_packages <- c("ggplot2", "readr", "grid")
#new_packages <- all_packages[!(all_packages %in% installed.packages()[,"Package"])]
#if(length(new_packages)) install.packages(new_packages)


# Load packages
library(ggplot2)
library(readr)
library(grid)

# Create custom theme for plots
customTheme = theme_bw(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=14, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(margin=margin(7,0,0,0), size=13),
        axis.title.y = element_text(size=14,  margin=margin(0,15,0,0)),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=11), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=13, margin=margin(-20,0,0,0))) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) 

# Set up the report
pdf(reportFile)
grid.newpage()
cover <- textGrob("TranscriptClean Report", gp=gpar(fontsize=30, col="black"))
grid.draw(cover)
data = read_delim(logFile, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)

# Plot 1: Size distribution of deletions
deletions = subset(data, ErrorType == "Deletion")
ggplot(deletions, aes(factor(Size))) + customTheme + geom_bar(stat="count", fill="navy") + xlab("Deletion length (bp)")  + ylab("Count") + theme(text= element_text(size=13)) + theme(axis.text.x = element_text(color = "black", size=12), axis.text.y = element_text(color = "black", size=12))
dev.off()

# Plot 1: Size distribution of insertions
insertions = subset(data, ErrorType == "Insertion")
ggplot(insertions, aes(factor(Size))) + customTheme + geom_bar(stat="count", fill="navy") + xlab("Insertion length (bp)")  + ylab("Count") + theme(text= element_text(size=13)) + theme(axis.text.x = element_text(color = "black", size=12), axis.text.y = element_text(color = "black", size=12))

# Plot 3: If noncanonical splice junction correction mode enabled, plot distribution of distance to nearest annotated junction






dev.off()






#function.name <- function(arg1, arg2, arg3=2, ...) {
#  newVar <- sin(arg1) + sin(arg2)  # do Some Useful Stuff
#  newVar / arg3   # return value 
#}
