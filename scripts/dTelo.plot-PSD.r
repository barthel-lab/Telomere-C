library(ggplot2)
library(gridExtra)
library(dplyr)
# Input
lig_bed <- read.delim("results/align/dTelo/BJ_capture.SES.bin.100000.window.bed", header = F, stringsAsFactors = F)
input_bed <- read.delim("results/align/dTelo/BJ_input.SES.bin.100000.window.bed", header = F, stringsAsFactors = F)
ITS_bed <- read.delim("data/chm13.v2.ITS.100K.bed", header = F, stringsAsFactors = F)

colnames(lig_bed) <- c("chr", "start", "end", "counts")
colnames(input_bed) <- c("chr", "start", "end", "counts")
colnames(ITS_bed) <- c("chr", "start", "end", "counts")


plist <- list()
qlist <- list()
# Plot for each chromosome
for(nchr in 1:22){
  sub_chr <- paste0("chr",nchr)
  lig_bed_sub <- subset(lig_bed, lig_bed$chr == sub_chr)
  input_bed_sub <- subset(input_bed, input_bed$chr == sub_chr)
  ITS_bed_sub <- subset(ITS_bed, ITS_bed$chr == sub_chr)
  delta_bed <- log(lig_bed_sub$counts / input_bed_sub$counts)

  # p_arm, 10M from end
  start_pos <- 1
  end_pos <- 200
  data <- data.frame("bin"=1:200, 
                     "delta"=delta_bed[start_pos:end_pos], 
                     "ITS"=ITS_bed_sub$counts[start_pos:end_pos])
  data <- na.omit(data)
  data <- data[!is.infinite(data$delta),]
  data["delta_percent"] <- data$delta/sum(data$delta)*100
  
  # Find periodicity by PSD. Spans can be none, 2, 2.5, 3.3, 5, or 10
  x.spec <- spectrum(data$delta_percent, spans=2, plot=FALSE)
  period <- 1/x.spec$freq[which.max(x.spec$spec)]
  max_period <- 200/period
  period_seq <- period*seq(1:max_period)
  
  # Calculate pearson correlation
  correlation <- cor(data$bin, data$delta_percent, method = 'pearson')
  peasonR <- paste("r = ", round(correlation,4))
  
  # Plot
  plist[[nchr]] <- ggplot(data, aes(x=bin, y=delta_percent)) + 
    geom_point() + 
    geom_smooth(method = "glm") +
    labs(title = paste0(sub_chr,"p"),
         x = "distance to telomere * 100K",
         y = "log(Telo-C/input) %") + 
    geom_vline(xintercept = which(data$ITS>0), linetype="dotted", 
               color = "red", size=1) +
    annotate("text", x=median(data$bin), y=max(data$delta_percent), label=peasonR) + 
    geom_vline(xintercept = period_seq, linetype="dotted", size=1)
  
  #ggsave(filename = paste0(sub_chr,"p",".png"), p_arm, width=4, height=4)
  
  # q_arm, 10M from end
  start_pos <- length(delta_bed) 
  end_pos <- length(delta_bed) - 199
  data <- data.frame("bin"=1:200, 
                     "delta"=delta_bed[start_pos:end_pos],
                     "ITS"=ITS_bed_sub$counts[start_pos:end_pos])
  data <- na.omit(data)
  data <- data[!is.infinite(data$delta),]
  data["delta_percent"] <- data$delta/sum(data$delta)*100
  
  # Find periodicity by PSD. Spans can be none, 2, 2.5, 3.3, 5, or 10
  x.spec <- spectrum(data$delta_percent, spans=2, plot=FALSE)
  period <- 1/x.spec$freq[which.max(x.spec$spec)]
  max_period <- 200/period
  period_seq <- period*seq(1:max_period)
  
  # Calculate pearson correlation
  correlation <- cor(data$bin, data$delta_percent, method = 'pearson')
  peasonR <- paste("r = ", round(correlation,4))
  
  # Plot q_arm
  qlist[[nchr]] <- ggplot(data, aes(x=bin, y=delta_percent)) + 
    geom_point() + 
    geom_smooth(method = "glm") +
    labs(title = paste0(sub_chr,"q"),
         x = "distance to telomere  * 100K",
         y = "log(Telo-C/input) %") +
    geom_vline(xintercept = which(data$ITS>0), linetype="dotted", 
               color = "red", size=1) +
    annotate("text", x=median(data$bin), y=max(data$delta_percent), label=peasonR)+  
    geom_vline(xintercept = period_seq, linetype="dotted", size=1)
  
  #  +
  #    annotate("text", x=median(data$bin), y=max(data$delta_percent), label=Rsq) +
  #    annotate("text", x=median(data$bin), y=max(data$delta_percent)*0.95, label=my_eq)
  #ggsave(filename = paste0(sub_chr,"q",".png"), q_arm, width=4, height=4)
}

# manually remove low mapping region and short arms
plist <- plist[c(-13:-15,-18,-21,-22)]

ggsave(filename = "all_p_arm.png", arrangeGrob(grobs = plist, ncol = 4),width = 15, height = 15) 
ggsave(filename = "all_q_arm.png", arrangeGrob(grobs = qlist, ncol = 5),width = 15, height = 15) 
