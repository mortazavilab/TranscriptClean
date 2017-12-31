# Dana Wyman
# 12/29/2017
# This script generates a PDF report based on a TranscriptClean run
# Designed for use with R version R/3.3.2

main <-function() {

    # Read input arguments
    args = commandArgs(trailingOnly = TRUE)
    logFile = args[1]
    prefix = args[2]

    customTheme = setupRun()
    reportFile = paste(prefix, "report.pdf", sep="_")


    # Set up the report
    pdf(reportFile)
    grid.newpage()
    cover <- textGrob("TranscriptClean Report", gp=gpar(fontsize=28, col="black"))
    grid.draw(cover)


    # Read in data from run
    data = suppressMessages(read_delim(logFile, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE))

    # Plot 1: Size distribution of deletions
    # Median and max values are labeled on the plot
    deletions = subset(data, ErrorType == "Deletion")
    maxCount = max(table(factor(deletions$Size)))
    medianD = median(deletions$Size)
    lab = getMedMaxLabel(deletions$Size)

    p1 = ggplot(deletions, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") + 
    	xlab("Deletion length (bp)") + ylab("Count") + customTheme +
	ggtitle("Size distribution of deletions in transcripts prior to correction\n") + 
	geom_vline(aes(xintercept=medianD), color="grey", linetype="dashed", size=0.75) + 
	annotate("text", x = lineLabelPos(length(lab), medianD, max(deletions$Size)), y = maxCount*0.75, label = lab, color = "black")
    print(p1)

    # Plot 2: Size distribution of insertions
    # Median and max values are labeled on the plot
    insertions = subset(data, ErrorType == "Insertion")
    maxCount = max(table(factor(insertions$Size)))
    medianI = median(insertions$Size)
    lab = getMedMaxLabel(insertions$Size)
    p2 = ggplot(insertions, aes(Size)) + geom_bar(stat="count", fill="navy") +
    	xlab("Insertion length (bp)") + ylab("Count") + customTheme + 
	ggtitle("Size distribution of insertions in transcripts prior to correction\n") +
        geom_vline(aes(xintercept=medianD), color="grey", linetype="dashed", size=0.75) +
        annotate("text", x = lineLabelPos(length(lab), medianI, max(insertions$Size)), y = maxCount*0.75, label = lab, color = "black")
    print(p2)


    # Plot 3: If noncanonical splice junction correction mode enabled, plot distribution of distance to nearest annotated junction
    # Median and max values are labeled on the plot
    ncSJs = subset(data, ErrorType == "NC_SJ_boundary")
    maxCount = max(table(factor(ncSJs$Size)))
    medianS = median(ncSJs$Size)
    lab = getMedMaxLabel(ncSJs$Size)
    print(summary(ncSJs$Size))

    p3 = ggplot(ncSJs, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") +
        xlab("Distance from annotated splice site (bp)") + ylab("Count") + customTheme +
        ggtitle("Distribution of distance between noncanonical splice sites and \ntheir nearest annotated splice site\n") +
        geom_vline(aes(xintercept=medianS), color="grey", linetype="dashed", size=0.75) +
        annotate("text", x = lineLabelPos(length(lab), medianS, 100), y = maxCount*0.75, label = lab, color = "black") + 
        coord_cartesian(xlim = c(-50, 50))
    print(p3)

    # Plot 4: 


    dev.off()

}

# ------------------Functions------------------------------------------

setupRun <- function() {

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
    # axis.text controls tick mark labels
    customTheme = suppressMessages(theme_bw(base_family = "Helvetica", base_size = 14) +
        theme(plot.margin = unit(c(2.5,1,1,1), "cm")) +
        theme(plot.title = element_text(lineheight=1, size= 13.5, margin=margin(-10,1,1,1))) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +

        theme(axis.title.x = element_text(color="black", size=14, vjust = -2, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(color="black", vjust=0.75, size=13),
        axis.title.y = element_text(color="black", size=14, margin=margin(0,10,0,0)),
        axis.text.y  = element_text(color="black", vjust=0.75, size=13))  +
        theme(legend.text = element_text(color="black", size = 11), legend.title = element_text(color="black", size=11), legend.key.size = unit(0.5, "cm")))

    return(customTheme)
}

# ------------------Utilities for plotting------------------------------

lineLabelPos <- function(labelLen, linePos, axisLen) {
    # A common need in my plots is to be able to attach a label to a line
    # marking the median. The text position must take into account the 
    # overall width of the window to look good. This function figures out where
    # to put the label.
  
   offset = axisLen/10.0
   return(linePos + offset + labelLen)
}

getMedMaxLabel <- function(v) {
    # Returns a print-ready label of the median and max of vector v
    medianV = round(median(v), 2)
    maxV = max(v)
    medianLabel = paste("Median =", medianV, sep = " ")
    maxLabel = paste("Max =", maxV, sep = " ")
    plotLabel = paste( medianLabel, maxLabel, sep = "\n")
    return(plotLabel)
}

main()

#function.name <- function(arg1, arg2, arg3=2, ...) {
#  newVar <- sin(arg1) + sin(arg2)  # do Some Useful Stuff
#  newVar / arg3   # return value 
#}
