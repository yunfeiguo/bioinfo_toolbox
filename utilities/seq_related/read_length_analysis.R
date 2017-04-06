require('ggplot2')
require('flux')
require('data.table')
require('ggthemes')
require('grid')
require('gridExtra')
#scaling = 4.8
scaling = 1
setwd("/Users/guoy28/Desktop")
rawdata = read.table("homo_sequences_lens.txt")
rawdata = rawdata * scaling;
#median/quantiles
summary(rawdata$V1)
x = rawdata[,1]
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
breaks = seq(min(data1$len),max(data1$len),length.out = 300)
groupIndex = findInterval(data$len,breaks, all.inside = TRUE, left.open = TRUE, rightmost.closed = TRUE)
groupMidpoint = (breaks[groupIndex] + breaks[groupIndex + 1]) / 2
data1[, groupMidpoint := groupMidpoint]
data2[, groupMidpoint := groupMidpoint]
data1[, c("bp.sum") := list(sum(bp)), by = groupMidpoint]
data2[, c("bin.count") := list(sum(count)), by = groupMidpoint]
data1$panel = "bp count at given length"
data2$panel = "read count"
#unify variable for plotting
data1$y = data1$bp.sum
data2$y = data2$bin.count


data <- rbindlist(list(data1, data2), fill = TRUE)

#combined plots
ggplot(data = data, aes(x=len, y=y)) + facet_grid(panel~., scales = "free") + 
  geom_line(data = data1) + 
geom_bar(data = data2, stat = "identity") +
  xlab("read length (bp)") +
  scale_x_continuous(breaks = seq(0,7000*scaling,length.out = 8)) +
  ylab("y") +
  ggtitle("genia low-level simulation read length distribution") +
  geom_vline(aes(xintercept = mean(x), colour = "mean")) +
  geom_vline(aes(xintercept = median(x), colour = "median")) +
  geom_vline(aes(xintercept = n50, colour = "n50")) +
  #use scale_* to manually change mapping between color and strings!
  scale_colour_manual("legend: ", values=c(mean='red', median='blue', n50='orange')) + theme_bw() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5), text = element_text(size=20))

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


