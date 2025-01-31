getwd()
setwd("/Users/oashna/Documents/Shira_lab/PhD/thesis/review/scc2.scc1.ratio")

#calculate distance
cen <- read.table("cen.bed", stringsAsFactors = F)
bed$mid <- 1/2 * (bed$start + bed$end)
cen$mid <- 1/2 * (cen$V2 + cen$V3)
head(cen)
dist <- numeric(0)
for (i in c(1:966)){
  temp <- which(cen[,1] == bed[i,1])
  dist <- c(dist, (bed[i,21]-cen[temp,4]))
}
bed$dist <- dist
head(bed)
plot <- ggplot(bed, aes(x=dist, y=normWt)) +
        geom_point() +
        ylim(-0.5,6) +
  xlim(-500000,500000) +
  stat_summary_bin(fun='mean', bins=40,
                   color='red', linewidth=0.7, geom='line')
plot

plot2 <- ggplot(bed, aes(x=dist, y=normDD)) +
  geom_point(color = "red") +
  ylim(-0.5,6) +
  xlim(-500000,500000) +
  stat_summary_bin(fun='mean', bins=40,
                   color='blue', linewidth=0.7, geom='line')
plot2

combined_plot <- ggarrange(plot,
                           plot2,
                           nrow = 1,
                           ncol = 2)
combined_plot
