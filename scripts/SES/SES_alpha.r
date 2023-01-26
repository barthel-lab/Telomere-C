# Reference: Normalization, bias correction, and peak calling for ChIP-seq. 
# 2012 Stat Appl Genet Mol Biol
library("dplyr")
arg <- commandArgs(trailingOnly = T)

# Input
# $1 for ip_bed, $2 for input_bed
ip_bed <- read.delim(arg[1], header = F, stringsAsFactors = F)
input_bed <- read.delim(arg[2], header = F, stringsAsFactors = F)

print("Sorting data...")

# Put read counts side by side, then ascending order by ip reads counts
counts_sort <- data.frame(ip_bed$V4, input_bed$V4) %>% arrange(ip_bed.V4)

# Calculate cumulative sum, presented in percentage.
counts_sort["ip_percentCumSum"] <- cumsum(counts_sort$ip_bed.V4)/sum(counts_sort$ip_bed.V4)
counts_sort["input_percentCumSum"] <- cumsum(counts_sort$input_bed.V4)/sum(counts_sort$input_bed.V4)

# Calculate delta
counts_sort["delta"] <- abs(counts_sort$input_percentCumSum - counts_sort$ip_percentCumSum)

# find k-th bin, then calculate scaling factor:alpha
k <- which.max(counts_sort$delta)
ip_scaleFactor <- counts_sort$ip_percentCumSum[k]
input_scaleFactor <- counts_sort$input_percentCumSum[k]
alpha <- ip_scaleFactor/input_scaleFactor

# Output alpha
print(paste("The scaling factor is ", alpha))
write(alpha, file = "scale_factor.txt")

# Output SES plot and cut-off
print("Plotting SES plot and cut-off")
tiff(file="Signal_extraction_scaling.tiff")
plot(x=1:nrow(counts_sort)/nrow(counts_sort),y=counts_sort$ip_percentCumSum, col="blue",
     xlab = "% of bins",
     ylab = "% of tags")
par(new=TRUE)
plot(x=1:nrow(counts_sort)/nrow(counts_sort),y=counts_sort$input_percentCumSum,col="red",
     xlab = "% of bins",
     ylab = "% of tags")
abline(v=k/nrow(counts_sort), lty=2 , lwd=4,col="green")
dev.off()

