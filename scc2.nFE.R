getwd()
CBS=bed[,c(1:6)]      

#input bg file
wt2FE <- read.table("wt.bedGraph", stringsAsFactors = F)
dd2FE <- read.table("dd.bedGraph", stringsAsFactors = F)
dW2FE <- read.table("dWpl1.bedGraph", stringsAsFactors = F)

Et2FE <- read.table ("pds5_Et.bedGraph", stringsAsFactors = F)
IAA2FE <- read.table ("pds5_IAA.bedGraph", stringsAsFactors = F)

wtGlu2FE <- read.table ("wt_glu.bedGraph", stringsAsFactors = F)
GalEGlu2FE <- read.table("GalE_glu.bedGraph", stringsAsFactors = F)

#define a function: extract chrom number 
chrom_No <- function(df) {
  df <- df[df$V1 != "chr2micron",]
  temp <- strsplit(df[,1], "chr")
  temp2<- numeric(0)
  for(i in c(1:length(temp))){
    temp2<- c(temp2, as.numeric(as.roman(temp[[i]][2])))
  }
  df$chrom2 <- temp2
  assign('df',df, envir=.GlobalEnv) #or df <<- df
  head(df,5)
}

chrom_No(Et2FE)
Et2FE <- df
head(Et2FE)

#define a function: calculate max FE at CBS 
maxFE <- function(x, df){
  x <- numeric(0)
  temp <- numeric(0)
  for (i in c(1:966)){
    temp <- which(df[,5] == CBS[i,6] &
                    df[,2] >= (CBS[i,2]-50) &
                    df[,2] <= (CBS[i,3]-50))
    x <- c(x,max(df[temp,4]))
    assign("x",x,envir=.GlobalEnv)
  }
}

maxFE(x,Et2FE)
CBS$Et <- x
head(CBS)
median(CBS$wt)
median(CBS$dW)
median(CBS$dd)

#making df for boxplot
df1 <- CBS[,c(7:8)]
df1melt <- melt(df1,value.name="Scc2nFE")
df1melt$group <- "wt-dW-dd"
df1melt$variable <- factor(df1melt$variable,
                           levels=c("wt","dd"))

df2 <- CBS[,c(12,13)]
df2melt <- melt(df2, value.name="Scc2nFE")
df2melt$group <- "Gal-Eco1"
df2melt$variable <- factor(df2melt$variable,
                           levels=c("wtGlu","GalEGlu"))

df3 <- CBS[,c(10,11)]
df3melt <- melt(df3, value.name="Scc2nFE")
df3melt$group <- "Pds5-AID"
df3melt$variable <- factor(df3melt$variable,
                           levels=c("Et","IAA"))

df_scc2 <- rbind(df1melt,df2melt,df3melt)
df_scc2$group <- factor(df_scc2$group,
                        levels=c("wt-dW-dd","Gal-Eco1","Pds5-AID"))

#plot
library(ggsignif)
b1 <- ggplot(df1melt, aes(variable, Scc2nFE, color=variable)) +
  geom_boxplot(notch = T) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="bottom") 
b1
b2 <- ggplot(df2melt, aes(variable, Scc2nFE, color=variable)) +
  geom_boxplot(notch = T) +
  scale_y_continuous(limits = c(0,10)) +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Dark2") 
b2
b3 <- ggplot(df3melt, aes(variable, Scc2nFE, color=variable)) +
  geom_boxplot(notch = T) +
  scale_y_continuous(limits = c(0,10)) +
  theme(legend.position="bottom") +
  scale_color_brewer(palette="Dark2") 
b3
library(ggpubr)


combined_plot <- ggarrange(b1,
                           b2,
                           nrow = 1,
                           ncol = 2)
combined_plot
scc2nFEbox <- ggplot(df_scc2, aes(variable, Scc2nFE, fill=variable)) +
  geom_boxplot(notch = T) +
  facet_wrap(~group) +
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous(limits = c(1,5))

scc2nFEbox
