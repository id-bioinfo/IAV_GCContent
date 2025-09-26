library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(ggpubr)

#Fig. 4b
d1 <- read.table('annotated_sporadic_avian_genomic_gc_forR.tsv', header = T, sep = "\t")
# With transparency (right)
p1 <- ggplot(data=d1, aes(x=GC, group=Groups, fill=Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x="Genomic GC Content",y = "Density") # theme_ipsum() + 
p1

#Supplementary Fig. 4
d3 <- read.table('feature_8genes.tsv', header = T, sep = "\t")
#HA
gc <- ggplot(data=d3, aes(x=HA_gc, group=Groups, fill=Groups)) + ggtitle("HA d = 0.73") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(39, 45.5))
CpG <- ggplot(data=d3, aes(x=HA_CpG, group=Groups, fill=Groups)) + ggtitle("d = 1.82") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(7, 28))
GpC <- ggplot(data=d3, aes(x=HA_GpC, group=Groups, fill=Groups)) + ggtitle("d = 4.01") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(25, 52))
GpG <- ggplot(data=d3, aes(x=HA_GpG, group=Groups, fill=Groups)) + ggtitle("d = 3.01") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(38, 88), breaks = seq(38, 88, by = 10))#scale_x_continuous(limits = c(50, 80))
CpC <- ggplot(data=d3, aes(x=HA_CpC, group=Groups, fill=Groups)) + ggtitle("d = 5.39") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1, colour = "#FF0000")) + 
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(23, 52))
HA_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
HA_plot 


#NA
gc <- ggplot(data=d3, aes(x=NA_gc, group=Groups, fill=Groups)) + ggtitle("NA d = 1.24") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(40, 47))
CpG <- ggplot(data=d3, aes(x=NA_CpG, group=Groups, fill=Groups)) + ggtitle("d = 1.97") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(7.5, 27.5))
GpC <- ggplot(data=d3, aes(x=NA_GpC, group=Groups, fill=Groups)) + ggtitle("d = 5.28") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(27.5, 52.5))
GpG <- ggplot(data=d3, aes(x=NA_GpG, group=Groups, fill=Groups)) + ggtitle("d = 6.34") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(38, 88), breaks = seq(38, 88, by = 10))#scale_x_continuous(limits = c(55, 95))
CpC <- ggplot(data=d3, aes(x=NA_CpC, group=Groups, fill=Groups)) + ggtitle("d = 4.05") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(25, 55))
NA_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
NA_plot

#NP
gc <- ggplot(data=d3, aes(x=NP_gc, group=Groups, fill=Groups)) + ggtitle("NP d = 1.58") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(44, 49))
CpG <- ggplot(data=d3, aes(x=NP_CpG, group=Groups, fill=Groups)) + ggtitle("d = 2.51") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(17.5, 32.5))
GpC <- ggplot(data=d3, aes(x=NP_GpC, group=Groups, fill=Groups)) + ggtitle("d = 2.63") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(40, 55))
GpG <- ggplot(data=d3, aes(x=NP_GpG, group=Groups, fill=Groups)) + ggtitle("d = 4.48") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(38, 88), breaks = seq(38, 88, by = 10))#scale_x_continuous(limits = c(62.5, 87.5))
CpC <- ggplot(data=d3, aes(x=NP_CpC, group=Groups, fill=Groups)) + ggtitle("d = 2.43") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(30, 52.5))
NP_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
NP_plot

#PA
gc <- ggplot(data=d3, aes(x=PA_gc, group=Groups, fill=Groups)) + ggtitle("PA d = 1.50") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(40.5, 45.5))
CpG <- ggplot(data=d3, aes(x=PA_CpG, group=Groups, fill=Groups)) + ggtitle("d = 7.43") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(15, 32.5))
GpC <- ggplot(data=d3, aes(x=PA_GpC, group=Groups, fill=Groups)) + ggtitle("d = 3.42") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(35, 52.5))
GpG <- ggplot(data=d3, aes(x=PA_GpG, group=Groups, fill=Groups)) + ggtitle("d = 3.18") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(38, 88), breaks = seq(38, 88, by = 10))#scale_x_continuous(limits = c(51, 69))
CpC <- ggplot(data=d3, aes(x=PA_CpC, group=Groups, fill=Groups)) + ggtitle("d = 4.99") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(36, 52.5))
PA_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
PA_plot


#PB1
gc <- ggplot(data=d3, aes(x=PB1_gc, group=Groups, fill=Groups)) + ggtitle("PB1 d = 1.74") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1, colour = "#FF0000")) + 
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(40.5, 45.5))
CpG <- ggplot(data=d3, aes(x=PB1_CpG, group=Groups, fill=Groups)) + ggtitle("d = 6.83") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(11, 30))
GpC <- ggplot(data=d3, aes(x=PB1_GpC, group=Groups, fill=Groups)) + ggtitle("d = 3.00") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(33.5, 46.5))
GpG <- ggplot(data=d3, aes(x=PB1_GpG, group=Groups, fill=Groups)) + ggtitle("d = 2.27") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(38, 88), breaks = seq(38, 88, by = 10))#scale_x_continuous(limits = c(52.5, 67.5))
CpC <- ggplot(data=d3, aes(x=PB1_CpC, group=Groups, fill=Groups)) + ggtitle("d = 3.00") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(35, 52.5))
PB1_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
PB1_plot

#PB2
gc <- ggplot(data=d3, aes(x=PB2_gc, group=Groups, fill=Groups)) + ggtitle("PB2 d = 1.33") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(41, 47))
CpG <- ggplot(data=d3, aes(x=PB2_CpG, group=Groups, fill=Groups)) + ggtitle("d = 2.28") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(12.5, 30))
GpC <- ggplot(data=d3, aes(x=PB2_GpC, group=Groups, fill=Groups)) + ggtitle("d = 3.36") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(32.5, 50))
GpG <- ggplot(data=d3, aes(x=PB2_GpG, group=Groups, fill=Groups)) + ggtitle("d = 7.56") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(40, 85), breaks = seq(40, 85, by = 9))#scale_x_continuous(limits = c(57.5, 80))
CpC <- ggplot(data=d3, aes(x=PB2_CpC, group=Groups, fill=Groups)) + ggtitle("d = 3.58") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(30, 50))
PB2_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
PB2_plot

#M1
gc <- ggplot(data=d3, aes(x=M1_gc, group=Groups, fill=Groups)) + ggtitle("M1 d = 1.43") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(46.5, 52.5))
CpG <- ggplot(data=d3, aes(x=M1_CpG, group=Groups, fill=Groups)) + ggtitle("d = 7.78") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1, colour = "#FF0000")) +
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(17.5, 40))
GpC <- ggplot(data=d3, aes(x=M1_GpC, group=Groups, fill=Groups)) + ggtitle("d = 3.01") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(56, 76))
GpG <- ggplot(data=d3, aes(x=M1_GpG, group=Groups, fill=Groups)) + ggtitle("d = 6.42") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) + 
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(38, 88), breaks = seq(38, 88, by = 10))#scale_x_continuous(limits = c(57.5, 90.5))
CpC <- ggplot(data=d3, aes(x=M1_CpC, group=Groups, fill=Groups)) + ggtitle("d = 2.75") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(37.5, 60))
M1_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
M1_plot


#NS1
gc <- ggplot(data=d3, aes(x=NS1_gc, group=Groups, fill=Groups)) + ggtitle("NS1 d = 1.11") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="GC Content (%)",y = "Density") + scale_x_continuous(limits = c(40, 52), breaks = seq(40, 52, by = 3))#scale_x_continuous(limits = c(41.5, 49.5))
CpG <- ggplot(data=d3, aes(x=NS1_CpG, group=Groups, fill=Groups)) + ggtitle("d = 3.58") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="CpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(8, 38), breaks = seq(8, 38, by = 6))#scale_x_continuous(limits = c(15, 37.5))
GpC <- ggplot(data=d3, aes(x=NS1_GpC, group=Groups, fill=Groups)) + ggtitle("d = 6.01") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1,colour = "#FF0000")) +
  labs(x="GpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(30, 75), breaks = seq(30, 75, by = 9))#scale_x_continuous(limits = c(35, 60))
GpG <- ggplot(data=d3, aes(x=NS1_GpG, group=Groups, fill=Groups)) + ggtitle("d = 12.17") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1, colour = "#FF0000")) +
  labs(x="GpG Frequency (бы)",y = "") + scale_x_continuous(limits = c(38, 88), breaks = seq(38, 88, by = 10))#scale_x_continuous(limits = c(35, 80))
CpC <- ggplot(data=d3, aes(x=NS1_CpC, group=Groups, fill=Groups)) + ggtitle("d = 3.52") + 
  geom_density(adjust=1.5, alpha=.4) + theme(legend.position="none", plot.title = element_text(hjust = 1)) +
  labs(x="CpC Frequency (бы)",y = "") + scale_x_continuous(limits = c(25, 60), breaks = seq(25, 60, by = 7))#scale_x_continuous(limits = c(30, 57.5))
NS1_plot <- ggarrange(gc,CpG,GpC,GpG,CpC, align = "h", nrow = 1, ncol = 5) #
NS1_plot
  
#labels = c("A","B", "c"), , widths = c(3, 1.5, 2)
figure <- ggarrange(HA_plot,NA_plot,NP_plot,PA_plot,PB1_plot,PB2_plot,M1_plot,NS1_plot,labels = c("A","B","c","D","E","F","G","H"),align = "h", nrow = 8, ncol = 1) #
figure
