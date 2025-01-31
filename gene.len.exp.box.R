setwd("/Users/oashna/Documents/Shira_lab/PhD/thesis/review/gene.len.expression/")

##compare gene length between Scc1-associated CEN-oriented and TEL-oriented genes
scc1Gene <- read.table("scc1P.gene.valid.bed", stringsAsFactors = F)
head(scc1Gene)

scc1Gene.len <- scc1Gene[,c(11,12)]
colnames(scc1Gene.len) <- c("orientation","length")

scc1Gene.exp <- scc1Gene[,c(11,14)]
colnames(scc1Gene.exp) <- c("orientation","Rpo21_FE")

## plot_geneLength
library(ggsignif)
box1 <- ggplot(scc1Gene.len, aes(orientation, length, fill=orientation)) +
  geom_boxplot(notch = T) +
  scale_fill_brewer(palette="Set1") +
  geom_signif(comparisons=list(c("inward","outward")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "black",
              size = 0.4 )
              
box1

lenIn <- scc1Gene.len[scc1Gene.len$orientation == "inward",2] ##741
lenOut <- scc1Gene.len[scc1Gene.len$orientation == "outward",2] ##739
wilcox.test(lenIn,lenOut,paired = F) ##p-value= 0.3207

## plot_geneExpression
box2 <- ggplot(scc1Gene.exp, aes(orientation, Rpo21_FE, fill=orientation)) +
  geom_boxplot(notch = T) +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous(trans='log2') + 
  geom_signif(comparisons=list(c("inward","outward")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "black",
              size = 0.4 )

box2

##compare gene length and RPO21 FE of all CEN- and TEL-oriented genes
###Rpo21 bg: pol2bg
###cen: cen2
###input all genes (about 6000)
gene <- read.table("Scer_bed_5cols.txt", stringsAsFactors = F)
head(gene)

###orientation
  ###arm position
arm <- list()
for (i in c(1:5997)){
  for (j in c(1:16)) {
    if (gene[i,1] == cen2[j,1] & gene[i,3] <= cen2[j,2]){
      arm <- append(arm, "left")
    } else if (gene[i,1] == cen2[j,1] & gene[i,2] >= cen2[j,3]) {
      arm <- append(arm, "right")
    } else if (gene[i,1] == cen2[j,1] & 
               gene[i,3] >= cen2[j,3] &
               gene[i,2] <= cen2[j,2]) {
      arm <- append(arm,"throughCEN")
    }
  }
}
gene$arm <- unlist(arm)
nrow(gene[gene$arm == "right",])
nrow(gene[gene$arm == "left",])

  ###gene orientation
direction <- list()
for (i in c(1:5997)) {
  if (gene[i,5] == "+" & gene[i,6] == "left" |
      gene[i,5] == "-" & gene[i,6] == "right"){
    direction <- append(direction, "inward")
  } else if (gene[i,5] == "-" & gene[i,6] == "left" |
             gene[i,5] == "+" & gene[i,6] == "right"){
    direction <- append(direction, "outward")
  } else {
    direction <- append(direction, "invalid")
  }
}
gene$orient <- unlist(direction)
nrow(gene[gene$orient == "inward",])
nrow(gene[gene$orient == "outward",])

  ###gene length
gene$len <- gene$V3 - gene$V2

  ###calculate max rpo21FE at each gene
rpo21FE <- numeric(0)
temp <- numeric(0)
for (i in c(1:5997)){
  temp <- which(pol2bg[,1] == gene[i,1] &
                  pol2bg[,2] >= (gene[i,2]-50) &
                  pol2bg[,2] <= (gene[i,3]-50))
  rpo21FE <- c(rpo21FE,max(pol2bg[temp,4]))
}
gene$rpo21FE <- rpo21FE
head(gene)
gene_valid <- gene[gene$rpo21FE != "-Inf",]
  ###plot
box3 <- ggplot(gene_valid, aes(orient, len, fill=orient)) +
  geom_boxplot(notch = T) +
  scale_fill_brewer(palette="Set1") +
  geom_signif(comparisons=list(c("inward","outward")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "black",
              size = 0.4 )

box3

box4 <- ggplot(gene_valid, aes(orient, rpo21FE, fill=orient)) +
  geom_boxplot(notch = T) +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous(trans='log2') +
  geom_signif(comparisons=list(c("inward","outward")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "black",
              size = 0.4 )

gene_in <- gene[gene$orient == "inward",]
nrow(gene_in) #2994
gene_out <- gene[gene$orient == "outward",]
nrow(gene_out) #3003
