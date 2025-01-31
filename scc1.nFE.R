wt=read.table("wt.bedGraph",stringsAsFactors = F)
dd=read.table("dd.bedGraph", stringsAsFactors = F)
wtGlu=read.table("wt_glu.bedGraph",stringsAsFactors = F)
galE=read.table("GalE_glu.bedGraph",stringsAsFactors = F)

#extract chrom number
chrom_No(dd)
dd <- df

#calculate maxFE
maxFE(x,galE)
CBS$galE1 <- x
head(CBS)

#make df for boxplot
df4 <- CBS[,c(14,15)]
df4melt <- melt(df4,value.name="Scc1nFE")
df4melt$variable <- factor(df4melt$variable,
                           levels=c("wt1","dd1"))
b4 <- ggplot(df4melt, aes(variable, Scc1nFE, color=variable)) +
  geom_boxplot(notch = T) +
  scale_y_continuous(limits = c(0,16)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="bottom") 
b4

df5 <- CBS[,c(16,17)]
df5melt <- melt(df5,value.name="Scc1nFE")
df5melt$variable <- factor(df5melt$variable,
                           levels=c("wtGlu1","galE1"))
b5 <- ggplot(df5melt, aes(variable, Scc1nFE, color=variable)) +
  geom_boxplot(notch = T) +
  scale_y_continuous(limits = c(0,16)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="bottom") 
b5

combined_plot <- ggarrange(b4,
                           b5,
                           nrow = 1,
                           ncol = 2)
combined_plot
