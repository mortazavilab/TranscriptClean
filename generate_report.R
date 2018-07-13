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
    pdf(reportFile, paper ='letter', width = 7, height = 7)
    grid.newpage()
    cover <- textGrob("TranscriptClean Report", gp=gpar(fontsize=28, col="black"))
    grid.draw(cover)

    # Read in data from run
    print("Reading log files............")
    data = read_delim(logFileTE, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE, na = "NA")
    transcripts = read_delim(logFileVerbose, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE, na = "NA")
    transcripts[,3:12] <- sapply(transcripts[,3:12], as.numeric)

    # Table 1
    print("Creating tables..............")
    primary = c("Primary mapping", nrow(subset(transcripts, Mapping == "primary")))
    multi = c("Non-primary mapping", nrow(subset(transcripts, Mapping == "non-primary"))) 
    unmap = c("Unmapped", nrow(subset(transcripts, Mapping == "unmapped")))
    total = c("Total", nrow(transcripts))
    t1 = rbind(primary, rbind(multi, rbind(unmap, total)))
    title_t1 <- textGrob("Transcripts in Input",gp=gpar(fontface="bold", fontsize=13), vjust = -6)
    t1 = tableGrob(t1, rows = NULL, cols = c("Transcript type", "Count"))
    gt1 = gTree(children=gList(t1, title_t1))
    

    # Table 2
    processedTranscripts = subset(transcripts, Mapping == "primary")
    del = c(sum(processedTranscripts$corrected_deletions), sum(processedTranscripts$uncorrected_deletions), sum(processedTranscripts$variant_deletions))
    ins = c(sum(processedTranscripts$corrected_insertions), sum(processedTranscripts$uncorrected_insertions), sum(processedTranscripts$variant_insertions))
    mis = c(sum(processedTranscripts$corrected_mismatches), NA, sum(processedTranscripts$variant_mismatches))
    ncsj = c(sum(processedTranscripts$corrected_NC_SJs), sum(processedTranscripts$uncorrected_NC_SJs), NA)
    categories = c("Deletions", "Insertions", "Mismatches", "Noncanonical jns")
    t2 = cbind.data.frame(categories, rbind.data.frame(del, rbind.data.frame(ins, rbind.data.frame(mis, ncsj))))
    t2$total = rowSums(t2[,2:4], na.rm = TRUE)
    t2$total[t2$total == 0] = NA
    t2$percent = round((t2[,2])*100.0/t2[,5], 2)
 
    title_t2 = textGrob("Summary of Processed Errors",gp=gpar(fontface="bold", fontsize=13), vjust = -6)
    t2 = tableGrob(t2, rows = NULL, cols = c("Error Type", "Corrected", "Not Correctable", "Variant", "Total", "Percent Corrected"))
    gt2 = gTree(children=gList(t2, title_t2))    

    
    # Table 3
    processedTranscripts$totD = rowSums(processedTranscripts[,3:5], na.rm = TRUE)
    processedTranscripts$totI = rowSums(processedTranscripts[,6:8], na.rm = TRUE)    
    processedTranscripts$totM = rowSums(processedTranscripts[,9:10], na.rm = TRUE)
    processedTranscripts$totNCSJ = rowSums(processedTranscripts[,11:12], na.rm = TRUE)
    del = c(nrow(subset(processedTranscripts, totD > 0)), nrow(subset(processedTranscripts, totD > 0 & (uncorrected_deletions > 0 | variant_deletions > 0)))) 
    ins = c(nrow(subset(processedTranscripts, totI > 0)), nrow(subset(processedTranscripts, totI > 0 & (uncorrected_insertions > 0 | variant_insertions > 0))))
    mis = c(nrow(subset(processedTranscripts, totM > 0)), nrow(subset(processedTranscripts, totM > 0 & variant_mismatches > 0)))
    ncsj = c(nrow(subset(processedTranscripts, totNCSJ > 0)), nrow(subset(processedTranscripts, totNCSJ > 0 & uncorrected_NC_SJs > 0))) 
    categories = c("Deletions", "Insertions", "Mismatches", "Noncanonical jns")
    t3 = cbind.data.frame(categories, rbind.data.frame(del, rbind.data.frame(ins, rbind.data.frame(mis, ncsj))))
    t3$change = round((t3[,2] - t3[,3])*100.0/t3[,2], 2)
    title_t3 = textGrob("Transcripts containing one or more of error/variant type",gp=gpar(fontface="bold", fontsize=13), vjust = -6)
    t3 = tableGrob(t3, rows = NULL, cols = c("Type", "Before TranscriptClean", "After TranscriptClean", "Percent Corrected"))
    gt3 = gTree(children=gList(t3, title_t3))

    # Drawing tables
    grid.arrange(gt1, gt2, gt3)#, layout_matrix = cbind(c(1,2,3),c(1,4,4)))


    # Plot 1: Size distribution of deletions
    # Median and max values are labeled on the plot
    print("Plot 1..................")
    deletions = subset(data, ErrorType == "Deletion")
    maxCount = max(table(factor(deletions$Size)))
    medianD = median(deletions$Size)
    lab = getMedMaxLabel(deletions$Size)
    q = quantile(deletions$Size, probs = 0.99)[1]

    p1 = ggplot(deletions, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") + 
    	xlab("Deletion length (bp)") + ylab("Count") + customTheme +
	ggtitle("Size distribution of deletions in transcripts prior to correction\n") +
        coord_cartesian(xlim = c(0,2*q))  +
	annotate("text", x = q, y = maxCount*0.75, label = lab, color = "black", size = 6)
    print(p1)
    
    # Plot 2: Size distribution of insertions
    # Median and max values are labeled on the plot
    print("Plot 2..................")
    insertions = subset(data, ErrorType == "Insertion")
    maxCount = max(table(factor(insertions$Size)))
    medianI = median(insertions$Size)
    lab = getMedMaxLabel(insertions$Size)
    q = quantile(insertions$Size, probs = 0.99)[1]

    p2 = ggplot(insertions, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") +
    	xlab("Insertion length (bp)") + ylab("Count") + customTheme + 
	ggtitle("Size distribution of insertions in transcripts prior to correction\n") +
        coord_cartesian(xlim = c(0,2*q)) +
        annotate("text", x = q, y = maxCount*0.75, label = lab, color = "black", size = 6)
    print(p2)
    
    # Plot 3: Size distribution of all indels
    # Median and max values are labeled on the plot
    print("Plot 3..................")
    indels = subset(data, ErrorType == "Deletion" | ErrorType == "Insertion")
    maxCount = max(table(factor(indels$Size)))
    lab = getMedMaxLabel(indels$Size)
    q = quantile(indels$Size, probs = 0.99)[1]
    lab = paste(lab, "\n", "99% Quantile = ", q, sep="")

    p3 = ggplot(indels, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") +
        xlab("Indel length (bp)") + ylab("Count") + customTheme +
        ggtitle("Size distribution of indels in transcripts prior to correction\n") +
        coord_cartesian(xlim = c(0,2*q)) +
        annotate("text", x = q, y = maxCount*0.75, label = lab, color = "black", size = 6)
    print(p3)
    
    # Plot 4: If noncanonical splice junction correction mode enabled, plot distribution of distance to nearest annotated junction
    # Median and max values are labeled on the plot
    ncSJs = subset(data, ErrorType == "NC_SJ_boundary")
    if (nrow(ncSJs) != 0) {
        print("Plot 4..................")
        ncSJs$Size = abs(ncSJs$Size)
        maxCount = max(table(factor(ncSJs$Size)))
        medianS = median(ncSJs$Size)
        lab = getMedMaxLabel(ncSJs$Size)
        # Find the mode
        mode = names(table(ncSJs$Size))[table(ncSJs$Size)==max(table(ncSJs$Size))]
        lab = paste(lab, "\n", "Mode = ", mode, sep = "")

        p4 = ggplot(ncSJs, aes(Size)) + geom_bar(stat="count", fill="dodgerblue4") +
            xlab("Distance from annotated splice site (bp)") + ylab("Count") + customTheme +
            ggtitle("Distribution of distance between noncanonical splice sites and \ntheir nearest annotated splice site\n") +
            geom_vline(aes(xintercept=medianS), color="grey", linetype="dashed", size=0.75) +
            annotate("text", lineLabelPos(1, medianS, medianS*2), y = maxCount*0.75, label = lab, color = "black", size = 6) + 
            coord_cartesian(xlim = c(0, medianS*2))
        print(p4)
    }
     
    # Plot 5: Overview of corrections made to insertions, deletions, mismatches, and noncanonical splice sites
    if (nrow(subset(data, Corrected == "Corrected")) > 0) {
        print("Plot 5..................")
        data_p5 = data
        data_p5[data_p5$ErrorType == "NC_SJ_boundary", "ErrorType"] = "NC_SJ"
        data_p5$Category = rep("Corrected", nrow(data_p5))
        catOrder = c("Corrected")
        plotcolors = c("skyblue")
        if (nrow(subset(data_p5, ReasonNotCorrected == "VariantMatch")) > 0) {
            data_p5$Category[data_p5$ReasonNotCorrected == "VariantMatch"] <- "Variant"
            catOrder = c(catOrder, "Variant")
            plotcolors = c(plotcolors, "navy")
        }
        if (nrow(subset(data_p5, ReasonNotCorrected == "TooLarge")) > 0) {
            data_p5$Category[data_p5$ReasonNotCorrected == "TooLarge"] <- "Uncorrected (Too large)"
            catOrder = c(catOrder, "Uncorrected (Too large)")
            plotcolors = c(plotcolors, "red")
        }
        if (nrow(subset(data_p5, ReasonNotCorrected == "TooFarFromAnnotJn")) > 0) {
            data_p5$Category[data_p5$ReasonNotCorrected == "TooFarFromAnnotJn"] <- "Uncorrected (Too far from annotated junction)"
            catOrder = c(catOrder, "Uncorrected (Too far from annotated junction)")
            plotcolors = c(plotcolors, "orange")   
        }
        if (nrow(subset(data_p5, ReasonNotCorrected == "Other")) > 0) {
            data_p5$Category[data_p5$ReasonNotCorrected == "Other"] <- "Uncorrected (Exon smaller than correction dist)"
            catOrder = c(catOrder, "Uncorrected (Exon smaller than correction dist)")
            plotcolors = c(plotcolors, "yellow")
        }

         p5 = ggplot(data_p5, aes(x=ErrorType, fill=factor(Category, levels=catOrder))) + geom_bar(position = "dodge") +
             xlab("") + ylab("Count") + customTheme + scale_fill_manual("",values = plotcolors) +
            ggtitle("Overview of corrections made to insertions, deletions, \nmismatches, and noncanonical splice sites") +
            theme(axis.text.x = element_text(angle = 0), legend.position = "bottom", legend.direction = "vertical")
        print(p5)
    }    
    
    # Plot 6: Percentage of transcripts containing error of a given type before and after TranscriptClean 
    if (nrow(subset(data, Corrected == "Corrected")) > 0) {
        print("Plot 6..................")
        totalTranscripts = nrow(processedTranscripts)

        transcriptsDBefore = length(unique(subset(data, ErrorType == "Deletion")$TranscriptID)) 
        transcriptsDAfter = length(unique(subset(data, ErrorType == "Deletion" & Corrected == "Uncorrected")$TranscriptID))

        transcriptsIBefore = length(unique(subset(data, ErrorType == "Insertion")$TranscriptID))
        transcriptsIAfter = length(unique(subset(data, ErrorType == "Insertion" & Corrected == "Uncorrected")$TranscriptID))

        transcriptsMBefore = length(unique(subset(data, ErrorType == "Mismatch")$TranscriptID))
        transcriptsMAfter = length(unique(subset(data, ErrorType == "Mismatch" & Corrected == "Uncorrected")$TranscriptID))

        transcriptsSJBefore = length(unique(subset(data, ErrorType == "NC_SJ_boundary")$TranscriptID))
        transcriptsSJAfter = length(unique(subset(data, ErrorType == "NC_SJ_boundary" & Corrected == "Uncorrected")$TranscriptID))

        data_p6 = data.frame(ErrorType=c("Deletion", "Insertion", "Mismatch", "NC_SJ"), Before=c(transcriptsDBefore, transcriptsIBefore, transcriptsMBefore, transcriptsSJBefore), After=c(transcriptsDAfter, transcriptsIAfter, transcriptsMAfter, transcriptsSJAfter))
        data_p6 = suppressMessages(melt(data_p6))
        data_p6$percent = vapply(as.numeric(data_p6$value), percent, numeric(1), totalTranscripts)
        data_p6$percent = paste(data_p6$percent,"%",sep="")

        p6 = ggplot(data=data_p6, aes(x=ErrorType, y=value, fill = variable)) + geom_bar(stat="identity",position="dodge") +
            xlab("Error Type") + ylab("Number of transcripts containing at least one error") + customTheme + scale_fill_manual("", values = c("skyblue", "navy")) + 
            ggtitle("Transcripts containing a given error type before and after \nTranscriptClean correction") + 
            geom_text(data=data_p6, aes(label=percent), position=position_dodge(width=0.9), vjust=-0.25)
        print(p6)
    }
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
    library(gridExtra)
    library(reshape2)
    library(gtable)

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
