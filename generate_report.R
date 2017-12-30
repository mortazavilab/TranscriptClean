# Dana Wyman
# 12/29/2017
# This script generates a PDF report based on a TranscriptClean run
# Designed for use with R version _

# Read input arguments
args = commandArgs(trailingOnly = TRUE)
logFile = args[1]
prefox = args[2]
reportFile = paste(prefix, "report.pdf", sep="_")


# Check for packages and install if not found
all_packages <- c("ggplot2", "readr")
new_packages <- all_packages[!(all_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)


# Load packages
library(ggplot2)
library(readr)

