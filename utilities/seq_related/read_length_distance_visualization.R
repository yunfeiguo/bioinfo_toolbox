#!/usr/bin/env Rscript
require('ggplot2')
require('flux')
require('data.table')
require('ggthemes')
require('grid')
require('reshape')
require('gridExtra')

get_ticks <- function(x) {
  #get nice looking ticks for a sequence of integers
  #find appropriate intervals
  min = floor(min(x) * 0.8)
  max = ceiling(max(x) * 1.2)
  interval <- floor((max - min) / 10)
  interval <- interval + 10 - (interval %% 10)
  return(seq.int(min, max, by = interval))
}

args <- commandArgs(trailingOnly = TRUE)

scaling = 1
input = args[1]
title <- basename(input)
print(title)
rawdata = data.table(read.table(input, fill = TRUE, as.is = TRUE))
rawdata = rawdata[1:(which(rawdata$V1 == 'Total')[1]-1),]
rawdata = rawdata[,.(V1=as.integer(V1), V2 = as.integer(V2))]
rawdata = untable(data.table(rawdata[,V1]), num=rawdata$V2)
rawdata = rawdata * scaling;
#median/quantiles
summary(rawdata)
x = rawdata
x = sort(x, decreasing = TRUE)
#N50
n50=x[cumsum(x) >= 0.5*sum(x)][1]

#process data for cumulative bp count plot
n = length(x)
y = rep(0, length(unique(x)))
count = rep(0, length(unique(x)))
x = sort(x, decreasing = FALSE)
uniqueIndex = 0
for (i in 1:n) {
  if (i == 1 || x[i] != x[i-1]) {
    uniqueIndex = uniqueIndex + 1;
  }
  y[uniqueIndex] = y[uniqueIndex] + x[i]
  count[uniqueIndex] = count[uniqueIndex] + 1;
}

data1 = data.table(len=sort(unique(x), decreasing = FALSE),bp=y)
data2 = data.table(len=sort(unique(x), decreasing = FALSE), count)
data3 = data.table(len=sort(unique(x), decreasing = FALSE), cdf = cumsum(count)/sum(count))
breaks = seq(min(data1$len),max(data1$len),length.out = 300)
groupIndex = findInterval(data1$len,breaks, all.inside = TRUE, left.open = TRUE, rightmost.closed = TRUE)
groupMidpoint = (breaks[groupIndex] + breaks[groupIndex + 1]) / 2
data1[, groupMidpoint := groupMidpoint]
data2[, groupMidpoint := groupMidpoint]
data3[, groupMidpoint := groupMidpoint]
data1[, c("bp.sum") := list(sum(bp)), by = groupMidpoint]
data2[, c("bin.count") := list(sum(count)), by = groupMidpoint]
data3[, c("cdf") := list(median(cdf)), by = groupMidpoint]
data1$panel = "bp count at given length"
data2$panel = "read count"
data3$panel <- "CDF"
#unify variable for plotting
data1$y = data1$bp.sum
data2$y = data2$bin.count
data3$y <- data3$cdf


data <- rbindlist(list(data1, data2, data3), fill = TRUE)

pdf(paste(input,".pdf",sep=""), width = 12, height = 10)
#combined plots
ggplot(data = data, aes(x=len, y=y)) + facet_grid(panel~., scales = "free") + 
  geom_line(data = data1) + 
geom_bar(data = data2, stat = "identity") +
  xlab("read length (bp)") +
  ylab("") + 
  scale_x_continuous(breaks = get_ticks(x)) +
  ggtitle(title) +
  geom_vline(aes(xintercept = mean(x), colour = "mean")) +
  geom_vline(aes(xintercept = median(x), colour = "median")) +
  geom_vline(aes(xintercept = n50, colour = "n50")) +
  #use scale_* to manually change mapping between color and strings!
  scale_colour_manual("legend: ", values=c(mean='red', median='blue', n50='orange')) + theme_bw() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5), text = element_text(size=20)) +
  geom_line(data = data3)

dev.off()


#separate plots
#histogram
if (FALSE) {
histPlot <- ggplot(rawdata, aes(V1))+geom_histogram(binwidth = 10) + ggtitle("genia low-level simulation read length distribution") +
  xlab("read length (bp)") + scale_x_continuous(breaks = seq(0,7000*scaling,length.out = 8)) +
  geom_vline(aes(xintercept = mean(x), colour = "mean")) +
  geom_vline(aes(xintercept = median(x), colour = "median")) +
  geom_vline(aes(xintercept = n50, colour = "n50")) +
  #use scale_* to manually change mapping between color and strings!
  scale_colour_manual("legend: ", values=c(mean='red', median='blue', n50='orange')) + theme_bw() +
  theme(legend.position = "top")

# throughput plot
throughputPlot <- ggplot(data1, aes(groupMidpoint, bp.sum)) + geom_line(aes(x=groupMidpoint, y=bp.sum)) +
  ggtitle("bp count vs read length for genia low-level simulation data") +
  xlab("read length") + scale_x_continuous(breaks = seq(0,7000*scaling,length.out = 8)) +
  ylab("bp count at given length") +
  geom_vline(aes(xintercept = n50, colour="n50")) +
               scale_colour_manual("legend: ", values = c(n50 = "orange")) +
  theme_bw() + theme(legend.position = "top")

#validation of N50 calculation
auc(data1[data1$groupMidpoint<n50,groupMidpoint], data1[data1$groupMidpoint<n50, bp.sum])/auc(data1$groupMidpoint, data1$bp.sum)


#combine two plots
grid.arrange(histPlot, throughputPlot, nrow = 2, ncol = 1)
}


