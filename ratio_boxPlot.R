getwd()


#input: bg file, bed file
bed <- read.table("36_sm1000_nocen.bed", stringsAsFactors = F)
wt1 <- read.table("scc2.scc1.wt.ratio.bedGraph", stringsAsFactors = F)
dd <- read.table("scc2.scc1.dd.ratio.bedGraph", stringsAsFactors = F)

#extract chrom name as number
  #for bed file
head(bed,10)
temp <- as.numeric(as.roman(bed$V1))
colnames(bed) <- c("chrom","start","end","width","peak")
bed$chrom2 <- temp

  #for bg files
    #wt
head(wt1,10)
wt12 <- wt1[wt1$V1 != "chr2micron",]
head(wt12,10)

temp<- strsplit(wt12[,1], "chr")
temp2<- numeric(0)
for(i in c(1:length(temp)){
  temp2<- c(temp2, as.numeric(as.roman(temp[[i]][2])))
}
table(temp2)
wt12$chr2 <- temp2
head(wt12,3)

    #dd
head(dd,10)
dd2 <- dd[dd$V1 != "chr2micron",]
head(dd2,10)

temp<- strsplit(dd2[,1], "chr")
temp2<- numeric(0)
for(i in c(1:length(temp)){
  temp2<- c(temp2, as.numeric(as.roman(temp[[i]][2])))
}
table(temp2)
dd2$chr2 <- temp2
head(dd2,3)

#calculate max ratio/peak region
  #for wt
wt1_ratio <- numeric(0)
temp <- numeric(0)
for (i in c(1:nrow(bed))){
  temp <- which(wt12[,5] == bed[i,6] &
                wt12[,2] >= (bed[i,2]-50) &
                wt12[,2] <= (bed[i,3]-50))
  wt1_ratio <- c(wt1_ratio,max(wt12[temp,4]))
}
bed$wt1_ratio <- wt1_ratio
head(bed)
  
  #for dd
dd_ratio <- numeric(0)
temp <- numeric(0)
for (i in c(1:nrow(bed)){
  temp <- which(dd2[,5] == bed[i,6] &
                  dd2[,2] >= (bed[i,2]-50) &
                  dd2[,2] <= (bed[i,3]-50))
  dd_ratio <- c(dd_ratio,max(dd2[temp,4]))
}
bed$dd_ratio <- dd_ratio
head(bed)

#take 950 random bin from non-binding regions
wt_nonbinding <- read.table("wt_ratio.nonbinding.nocen.bedGraph", stringsAsFactors = F)
dd_nonbinding <- read.table("dd_ratio.nonbinding.nocen.bedGraph", stringsAsFactors = F)

wt_nonbinding_rd <- wt_nonbinding[sample(nrow(wt_nonbinding),966),]
head(wt_nonbinding_rd,5)

dd_nonbinding_rd <- dd_nonbinding[sample(nrow(dd_nonbinding), 966),]
head(dd_nonbinding_rd,5)

#make dataframe for boxplot
bed$wt_nonbd <- wt_nonbinding_rd$V4
bed$dd_nonbd <- dd_nonbinding_rd$V4
head(bed,5)
bed_boxP <- bed[,c(7,8,9,10)]
head(bed_boxP,5)
bed_boxP1 <- bed[,c(7,8)]

install.packages(c("MASS","reshape2","reshape"))
library(MASS)
library(reshape)
library(reshape2)
bed2 <- bed[,c(7:10)]
bed2m <- melt(bed2, value.name = "log2_ratio")

bed_boxP2 <- melt(bed_boxP, na.rm=F, value.name = "ratio")
head(bed_boxP2,5)

bed_boxP1_melt <- melt(bed_boxP1, na.rm = F, value.name = "ratio")
#boxplot
library(ggplot2)
bed_boxP2$variable <- factor(bed_boxP2$variable, 
                             levels=c("wt_nonbd","dd_nonbd","wt_ratio","dd_ratio"))
bed2m$variable <- factor(bed2m$variable,
                         levels=c("wt1_ratio","dd_ratio","wt_ratio","galE_ratio"))
box <- ggplot(bed2m,aes(variable, log2_ratio)) +
  geom_boxplot(notch = T)
box

#normalization
median(bed2$wt1_ratio)

bed$rWtGal <- 2**bed$wt_ratio 
bed$rGalE <- 2**bed$galE_ratio
bed$rWt <- 2**bed$wt1_ratio
bed$rDD <- 2**bed$dd_ratio
median(bed$rWt)

head(bed)
log2((2**(median(bed$wt_ratio)))*(1/2**(median(bed$wt_ratio))))

bed$normWtGal <- log2((2**bed$wt_ratio) * (1/median(2**bed$wt_ratio)))
bed$normGalE <- log2(bed$rGalE * (1/median(2**bed$wt_ratio)))
bed$normDD <- log2(bed$rDD * (1/0.7383672))
bed$normWt <- log2((2**bed$wt1_ratio) * (1/median(2**bed$wt1_ratio)))
median(bed$normWtGal)
median(bed$normGalE)
median(bed$normWt)
median(bed$normDD)
2**0.3724064
2**1.342164

nbed1 <- bed[,c(13:14)]
colnames(nbed1) <- c("treatment","control")
nbed1m <- melt(nbed1, value.name="norm_log2_ratio")
nbed1m$variable <- factor(nbed1m$variable,
                        levels=c("control","treatment"))
                         #,"normWtGal","normGalE"))
nbed1m$group <- "wt-dd"

nbed2 <- bed[,c(11:12)]
colnames(nbed2) <- c("control","treatment")
nbed2m <- melt(nbed2, value.name="norm_log2_ratio")
nbed2m$variable <- factor(nbed2m$variable,
                          levels=c("control","treatment"))
nbed2m$group <- "Gal-Eco1"

nbed3 <- rbind(nbed1m,nbed2m)
nbed3$group <- factor(nbed3$group, 
                      levels = c("wt-dd","Gal-Eco1"))
box3 <- ggplot(nbed3, aes(variable, norm_log2_ratio, fill=variable)) +
  geom_boxplot(notch = T) +
  facet_wrap(~group) +
  scale_fill_brewer(palette="Dark2")
box3
box2 <- ggplot(nbedm,aes(variable, norm_log2_ratio, fill = variable)) +
  geom_boxplot(notch = T)
box2
box <- ggplot(bed_boxP2,aes(variable, ratio, fill=variable)) +
  geom_violin(alpha = 0.5)
box

box2 <- ggplot(bed_boxP2,aes(variable, ratio)) +
  geom_boxplot(notch = T)
box2

#Scc2 peaks bed file
bed2 <- read.table("Scc2-PK_wt_rep1.bed", stringsAsFactors = F)
head(bed2,10)
temp <- strsplit(bed2[,1], "chr")
length(temp)
temp2<- numeric(0)
for(i in c(1:length(temp)){
  temp2<- c(temp2, as.numeric(as.roman(temp[[i]][2])))
}
table(temp2)
bed2$chr2 <- temp2

#calculate max ration/Scc2 peak
  #wt
wt_ratio2 <- numeric(0)
temp <- numeric(0)
for (i in c(1:nrow(bed2)){
  temp <- which(wt2[,5] == bed2[i,6] &
                  wt2[,2] >= (bed2[i,2]-50) &
                  wt2[,2] <= (bed2[i,3]-50))
  wt_ratio2 <- c(wt_ratio2,max(wt2[temp,4]))
}
bed2$wt_ratio <- wt_ratio2
head(bed2,10)

  #dd
dd_ratio2 <- numeric(0)
temp <- numeric(0)
for (i in c(1:nrow(bed2)){
  temp <- which(dd2[,5] == bed2[i,6] &
                  dd2[,2] >= (bed2[i,2]-50) &
                  dd2[,2] <= (bed2[i,3]-50))
  dd_ratio2 <- c(dd_ratio2,max(dd2[temp,4]))
}
bed2$dd_ratio <- dd_ratio2
head(bed2,10)
colnames(bed2) <- c("chrom","start","end","name","peak","chr2","wt_ratio2","dd_ratio2")
bed2_box <- bed2[,c(7,8)]
bed2_box_melt <- melt(bed2_box, value.name="ratio")
head(bed2_box_melt)

bed_boxP3 <- rbind(bed_boxP2,bed2_box_melt)

#plot
box <- ggplot(bed_boxP3,aes(variable, ratio)) +
  geom_violin(alpha = 0.5)
box

box2 <- ggplot(bed_boxP3,aes(variable, ratio)) +
  geom_boxplot(notch = T)
box2

box1 <- ggplot(bed_boxP1_melt, aes(variable, ratio, fill=variable)) +
  geom_boxplot(notch = T)
box1 
