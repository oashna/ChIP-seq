getwd()
setwd("/Users/oashna/OneDrive - The University of Tokyo/MicroC/Scer.genes/fig1D.norm2/")
peak <- read.table("2020_004B_36.nocen.bed", stringsAsFactors = F)
ratio <- read.table("230607-PCE-sm1000_54I1.norm2.ratio.100.bedGraph", stringsAsFactors = F)
cpm <- read.table("54I1.100.bedGraph",stringsAsFactors = F)
wtratio <- read.table("220525-PCE-sm1000_80eIP.ratio.100.bedGraph", stringsAsFactors = F)
ratio36 <- read.table("2020_004B_36.ratio.100.bedGraph", stringsAsFactors = F)
ratio38 <- read.table("2020_004B_38.ratio.100.bedGraph", stringsAsFactors = F)

#ratio and cpm @peak
peak <- peak[,c(1:3)]

c <- numeric(0)
temp <- numeric(0)
for (i in c(1:nrow(peak2))) {
  temp <- which(pol2bg[,1] == peak2[i,1] & pol2bg[,2] >= (peak2[i,2] - 50) & pol2bg[,2] <= (peak2[i,3] -50))
  c <- c(c, max(pol2bg[temp,4]))
}
peak2$pol2bg <- c

#calculate distance
cen2 <- read.table("cen.bed", stringsAsFactors = F)
peak$mid <- 1/2 * (peak$V2 + peak$V3)
cen2$mid <- 1/2 * (cen2$V2 + cen2$V3)
dist <- numeric(0)
for (i in c(1:nrow(peak))){
  temp <- which(cen2[,1] == peak[i,1])
  dist <- c(dist, abs(peak[i,6]-cen2[temp,4]))
}
peak$dist2 <- dist
library("ggplot2")
peak2 <- peak[peak$ratio < 50,]
plot <- ggplot(peak2, aes(x=dist, y=wtratio)) +
  geom_point(aes(col=quant54)) +
  ylim(0,6) +
  stat_summary_bin(fun='mean', bins=20,
                   color='red', linewidth=0.7, geom='line') 
plot
peak2$quant54 <- quant_groups(peak2[,4],groups=4) 
peak2$dist2quant <- quant_groups(peak2[,11],groups=10)
plot <- ggplot(peak2, aes(x=ratio38, y=ratio)) +
  geom_point(aes(col=dist),size=0.8) +
  scale_color_viridis_c(option="inferno") 
library(ggpubr)
plot <- ggscatter(peak2, x="ratio36", y="wtratio",color ="dist2quant",
          add="reg.line", conf.int=T, cor.coef=T, 
          cor.method="pearson",cor.coef.size = 8) 
  scale_color_viridis_c(option="turbo") 
plot

peak2.sub1 <- peak2[,c(8,4,12)]
peak2.sub1.melt <- melt(peak2.sub1,
                       id.vars = "dist2quant", measure.vars = c("wtratio","ratio"))
peak2.sub1.melt$variable <- factor(peak2.sub1.melt$variable,levels = c("wtratio","ratio"))
box <- ggplot(peak2.sub1.melt, aes(x=dist2quant,y=value,fill=variable)) +
  geom_boxplot(notch = T, outlier.shape = NA)
box

peak2.sub2 <- peak2[,c(9,10,12)]
peak2.sub2.melt <- melt(peak2.sub2,
                        id.vars = "dist2quant", measure.vars = c("ratio36","ratio38"))

box <- ggplot(peak2.sub2.melt, aes(x=dist2quant,y=value,fill=variable)) +
  geom_boxplot(notch = T, outlier.shape = NA) #+
box
library(ggpubr)
ggscatter(peak2, x="pol2bg", y="wtratio", 
          add="reg.line", conf.int=T, 
          cor.coef=T, cor.method="spearman", 
          xlab="pol2wt", ylab="Scc2wt")
library(dvmisc)
peak2$scc2dd.quant <- quant_groups(peak2$ratio,groups=4)
table(peak2$scc2dd.quant)
q1 <- peak2[peak2$scc2dd.quant == "[0.702,1.28]",]
q2 <- peak2[peak2$scc2dd.quant == "(1.28,1.58]",]
q3 <- peak2[peak2$scc2dd.quant == "(1.58,2.07]",]
q4 <- peak2[peak2$scc2dd.quant == "(2.07,6.42]",]

plot(q1$dist, q1$ratio, pch = 1, cex = 0.5, ylim = c(0,6),
     xlim = c(-500000,500000),
     xlab = "Distance", 
     ylab = "Scc2-ddel_maxCPM",
     main = "Scc2-ddel max cpm vs Distance")
points(q2$dist, q2$ratio, pch = 1, cex = 0.5, col="blue")
points(q3$dist, q3$ratio, pch = 1, cex = 0.5, col="green")
points(q4$dist, q4$ratio, pch = 1, cex = 0.5, col="red")

plot(density(q1$dist),xlim=c(-500000,500000),ylim=c(0,0.0000035))
lines(density(q2$dist), col="blue")
lines(density(q3$dist),col="green")
lines(density(q4$dist),col="red")

peak2$scc1FC <- log2(peak2$ratio38) - log2(peak2$ratio36)
peak2$scc2FC <- log2(peak2$ratio) - log2(peak2$wtratio)
ggscatter(peak2, x="ratio36", y="scc1FC",
          color = "scc2dd.quant", palette = c("black","blue","green","red"))
#          add="reg.line", conf.int=T, 
#          cor.coef=T, cor.method="spearman", 
          xlab="Scc1dd", ylab="Scc2dd")
group <- list()
for (i in c(1:941)){
  if (peak2[i,14] < 0 & peak2[i,15] < 0){
    group <- append(group,"bothDown")
  } else if (peak2[i,14] < 0 & peak2[i,15] > 0) {
    group <- append(group,"scc1Down")
  } else if (peak2[i,14] > 0 & peak2[i,15] > 0) {
    group <- append(group,"bothUp")
  } else if (peak2[i,14] > 0 & peak2[i,15] < 0) {
    group <- append(group, "scc1Up")
  }
}
peak2$group <- unlist(group)
table(peak2$group)
peak2$group <- factor(peak2$group, levels=c("bothDown","scc1Down","scc1Up","bothUp"))
peak$dist2 <- dist
library(ggpubr)
x <- ggscatter(peak2, x="scc2FC", y="scc1FC",size=0.8) +
               #color="group") + 
               #palette = c("black","blue","green","red")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)
ggpar(x,xlim=c(-2,2),ylim = c(-2,2))
x        
          add="reg.line", conf.int=T, 
          cor.coef=T, cor.method="spearman")
x <- ggscatter(peak2, x="dist", y="ratio",size=0.8,
               color="group",
               palette = c("black","blue","green","red"))
ggpar(x,xlim=c(-500000,500000)) #,ylim = c(0,2))
g1 <- peak2[peak2$group == "bothDown",]
g2 <- peak2[peak2$group == "scc1Down",]
g3 <- peak2[peak2$group == "scc1Up",]
g4 <- peak2[peak2$group == "bothUp",]
plot(density(g1$dist),xlim=c(-500000,500000),ylim=c(0,0.0000055))
lines(density(g2$dist), col="blue")
lines(density(g4$dist),col="red")
write.table(g4[,c(1:3)], file="bothUp.bed", sep = "\t",
            quote = F, row.names = F, col.names = F)
#color = "scc2dd.quant", 
#palette = c("black","blue","green","red"),
library(reshape2)
peak2.1 <- peak2[,c(8,4,17)]
peak2.1melt <- melt(peak2.1,id.vars = "scc2dd.quant",
                    measure.vars = c("wtratio","ratio"))
box <- ggplot(peak2.1melt, aes(x=scc2dd.quant,y=value,fill=variable)) +
  geom_boxplot(notch = T)
box

peak2.2 <- peak2[,c(12,13,17)]
peak2.2melt <- melt(peak2.2,id.vars = "scc2dd.quant",
                    measure.vars = c("ratio75","ratio74"))
box <- ggplot(peak2.2melt, aes(x=scc2dd.quant,y=value,fill=variable)) +
  geom_boxplot(notch = T)
box

write.table(peak2[peak2$scc2dd.quant == "[0.702,1.28]",c(1:3)],
            file = "scc1P.Q1.bed",
            sep = "\t",quote = F,col.names = F,row.names = F)
write.table(peak2[peak2$scc2dd.quant == "(1.28,1.58]",c(1:3)],
            file = "scc1P.Q2.bed",
            sep = "\t",quote = F,col.names = F,row.names = F)
write.table(peak2[peak2$scc2dd.quant == "(1.58,2.07]",c(1:3)],
            file = "scc1P.Q3.bed",
            sep = "\t",quote = F,col.names = F,row.names = F)
write.table(peak2[peak2$scc2dd.quant == "(2.07,6.42]",c(1:3)],
            file = "scc1P.Q4.bed",
            sep = "\t",quote = F,col.names = F,row.names = F)
