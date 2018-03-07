# Dana Wyman
# 12/29/2017
# This script generates a PDF report based on a TranscriptClean run
# Designed for use with R version R/3.3.2

main <-function() {

    options(scipen=10000)

    # Read input arguments
    args = commandArgs(trailingOnly = TRUE)
    prefix = args[1]

    logFileTE = paste(prefix, "_clean.TE.log", sep = "")
    logFileVerbose = paste(prefix, "_clean.log", sep = "")

    customTheme = setupRun()
    reportFile = paste(prefix, "report.pdf", sep="_")


    # Set up the report
    #pdf(reportFile, paper='USr')
    pdf(reportFile, paper ='letter', width = 7, height = 7)
    grid.newpage()
    cover <- textGrob("TranscriptClean Report", gp=gpar(fontsize=28, col="black"))
    grid.draw(cover)


    # Read in data from run
    data = suppressMessages(read_delim(logFileTE, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE))

    # Plot 1: Size distribution of deletions
    # Median and max values are labeled on the plot
    deletions = subset(data, ErrorType == "Deletion")
    maxCount = max(table(factor(deletions$Size)))
    medianD = median(deletions$Size)
    lab = getMedMaxLabel(deletions$Size)

    p1 = ggplot(deletions, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") + 
    	xlab("Deletion length (bp)") + ylab("Count") + customTheme +
	ggtitle("Size distribution of deletions in transcripts prior to correction\n") +
        coord_cartesian(xlim = c(0,10)) + scale_x_continuous(breaks=seq(1,10))  +
	annotate("text", x = 5, y = maxCount*0.75, label = lab, color = "black", size = 6)
    print(p1)
    
    # Plot 2: Size distribution of insertions
    # Median and max values are labeled on the plot
    insertions = subset(data, ErrorType == "Insertion")
    maxCount = max(table(factor(insertions$Size)))
    medianI = median(insertions$Size)
    lab = getMedMaxLabel(insertions$Size)
    p2 = ggplot(insertions, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") +
    	xlab("Insertion length (bp)") + ylab("Count") + customTheme + 
	ggtitle("Size distribution of insertions in transcripts prior to correction\n") +
        coord_cartesian(xlim = c(0,10)) + scale_x_continuous(breaks=seq(1,10))  +
        annotate("text", x = 5, y = maxCount*0.75, label = lab, color = "black", size = 6)
    print(p2)
    
    # Plot 3: Size distribution of all indels
    # Median and max values are labeled on the plot
    indels = subset(data, ErrorType == "Deletion" | ErrorType == "Insertion")
    maxCount = max(table(factor(indels$Size)))
    lab = getMedMaxLabel(indels$Size)
    p3 = ggplot(indels, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") +
        xlab("Indel length (bp)") + ylab("Count") + customTheme +
        ggtitle("Size distribution of indels in transcripts prior to correction\n") +
        coord_cartesian(xlim = c(0,10)) + scale_x_continuous(breaks=seq(1,10)) +
        annotate("text", x = 5, y = maxCount*0.75, label = lab, color = "black", size = 6)
    print(p3)
    
    # Plot 4: If noncanonical splice junction correction mode enabled, plot distribution of distance to nearest annotated junction
    # Median and max values are labeled on the plot
    ncSJs = subset(data, ErrorType == "NC_SJ_boundary")
    if (nrow(ncSJs) != 0) {
        maxCount = max(table(factor(ncSJs$Size)))
        medianS = median(ncSJs$Size)
        lab = getMedMaxLabel(ncSJs$Size)

        p4 = ggplot(ncSJs, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") +
            xlab("Distance from annotated splice site (bp)") + ylab("Count") + customTheme +
            ggtitle("Distribution of distance between noncanonical splice sites and \ntheir nearest annotated splice site\n") +
            geom_vline(aes(xintercept=medianS), color="grey", linetype="dashed", size=0.75) +
            annotate("text", x = lineLabelPos(length(lab), medianS, 100), y = maxCount*0.75, label = lab, color = "black") + 
            coord_cartesian(xlim = c(-50, 50))
        print(p4)
    }
    
    # Plot 5: Overview of corrections made to insertions, deletions, mismatches, and noncanonical splice sites
    if (nrow(subset(data, Corrected == "Corrected")) > 0) {
        data[data$Corrected == "True", "Corrected"] = "Corrected"
        data[data$Corrected == "False", "Corrected"] = "Uncorrected"
        data_p4 = within(data, ReasonNotCorrected <- paste('(',ReasonNotCorrected, ')', sep=''))
        data_p4 = within(data_p4, Category <- paste(ErrorType,Corrected,ReasonNotCorrected,sep=' '))
        data_p4$Category <- gsub(' \\(NA\\)', '', data_p4$Category)
        plotcolors = c("red1", "red4", "red3", "goldenrod1", "darkorange4", "darkorange", "springgreen3", "springgreen4", "skyblue", "navy")
        catOrder = c("Deletion Uncorrected (TooLarge)", "Deletion Uncorrected (VariantMatch)", "Deletion Corrected",
                     "Insertion Uncorrected (TooLarge)", "Insertion Uncorrected (VariantMatch)", "Insertion Corrected",
                     "Mismatch Uncorrected (VariantMatch)", "Mismatch Corrected",
                     "NC_SJ_boundary Uncorrected (TooFarFromAnnotJn)", "NC_SJ_boundary Corrected")   

         p5 = ggplot(data_p4, aes(x=ErrorType, fill=factor(Category, levels=catOrder))) + geom_bar() +
             xlab("") + ylab("Count") + customTheme + scale_fill_manual("",values = plotcolors) +
            ggtitle("Overview of corrections made to insertions, deletions, \nmismatches, and noncanonical splice sites") +
            theme(axis.text.x = element_text(angle = 45))
        print(p5)
    }    
    quit()
    # Plot 5: Percentage of transcripts containing error of a given type before and after TranscriptClean 
    cmd = paste("wc -l <", logFileVerbose, sep = " ") 
    totalTranscripts = as.numeric(system(cmd, intern = TRUE)) # get total transcript number from other log file because TE log only records errors
    print(totalTranscripts)

    transcriptsDBefore = length(unique(subset(data, ErrorType == "Deletion")$TranscriptID)) 
    transcriptsDAfter = length(unique(subset(data, ErrorType == "Deletion" & Corrected == "Uncorrected")$TranscriptID))

    transcriptsIBefore = length(unique(subset(data, ErrorType == "Insertion")$TranscriptID))
    transcriptsIAfter = length(unique(subset(data, ErrorType == "Insertion" & Corrected == "Uncorrected")$TranscriptID))

    transcriptsMBefore = length(unique(subset(data, ErrorType == "Mismatch")$TranscriptID))
    transcriptsMAfter = length(unique(subset(data, ErrorType == "Mismatch" & Corrected == "Uncorrected")$TranscriptID))

    transcriptsSJBefore = length(unique(subset(data, ErrorType == "NC_SJ_boundary")$TranscriptID))
    transcriptsSJAfter = length(unique(subset(data, ErrorType == "NC_SJ_boundary" & Corrected == "Uncorrected")$TranscriptID))

    data_p5 = data.frame(ErrorType=c("Deletion", "Insertion", "Mismatch", "NC_SJ"), Before=c(transcriptsDBefore, transcriptsIBefore, transcriptsMBefore, transcriptsSJBefore), After=c(transcriptsDAfter, transcriptsIAfter, transcriptsMAfter, transcriptsSJAfter))
    data_p5 = melt(data_p5)
    data_p5$percent = vapply(as.numeric(data_p5$value), percent, numeric(1), totalTranscripts)
    data_p5$percent = paste(data_p5$percent,"%",sep="")

    p5 = ggplot(data=data_p5, aes(x=ErrorType, y=value, fill = variable)) + geom_bar(stat="identity",position="dodge") +
        xlab("Error Type") + ylab("Number of transcripts containing at least one error") + customTheme + scale_fill_manual("", values = c("skyblue", "navy")) + 
        ggtitle("Transcripts containing a given error type before and after \nTranscriptClean correction") + 
        geom_text(data=data_p5, aes(label=percent), position=position_dodge(width=0.9), vjust=-0.25)
    print(p5)
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
    library(reshape2)

    # Create custom theme for plots
    # axis.text controls tick mark labels
    customTheme = suppressMessages(theme_bw(base_family = "Helvetica", base_size = 14) +
        theme(plot.margin = unit(c(2.5,1,1,1), "cm")) +
        theme(plot.title = element_text(lineheight=1, size= 13.5, margin=margin(-10,1,1,1))) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +

        theme(axis.title.x = element_text(color="black", size=18, vjust = -2, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(color="black", vjust=0.75, size=16),
        axis.title.y = element_text(color="black", size=18, margin=margin(0,10,0,0)),
        axis.text.y  = element_text(color="black", vjust=0.75, size=16))  +
        theme(legend.text = element_text(color="black", size = 14), legend.title = element_text(color="black", size=11), legend.key.size = unit(0.5, "cm")))

    return(customTheme)
}

# ------------------Utilities for plotting------------------------------

percent <- function(x, tot) {
    # Calculate the percentage of x relative to total, and round to two decimal places
    p = round(100*x/tot, 2)
    return(p)
}

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
    v = as.numeric(v)
    medianV = round(median(v), 3)
    meanV = round(mean(v), 3)
    maxV = max(v)
    medianLabel = paste("Median =", medianV, sep = " ")
    meanLabel = paste("Mean =", meanV, sep = " ")
    maxLabel = paste("Max =", maxV, sep = " ")
    plotLabel = paste( medianLabel, meanLabel, maxLabel, sep = "\n")
    return(plotLabel)
}

main()
