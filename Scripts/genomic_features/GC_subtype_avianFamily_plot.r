library(ggplot2)

#Fig. 4a
data1 = read.csv(file = 'GC_subtype_avianFamily_data/av_spo_genomicGC.csv')
p1 <- ggplot(data1, aes(x=Group, y=GenomicGC, fill=Label)) +
  geom_violin(width=1, size=0.4) + 
  scale_fill_manual(values=c("#DCDCDC","#778899")) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outliers = FALSE) +
  #geom_jitter(height = 0, size = 0.5,  width = 0.05, aes(colour = factor(Preidct))) +
  theme_classic() + ylim(42, 46)
p1

#
data2 = read.csv(file = 'GC_subtype_avianFamily_data/per_genomicGC_scaledTime.csv')
p2 <- ggplot(data2, aes(factor(Group), GenomicGC)) +
  geom_violin(width=1, size=0.2) + theme_classic() +
  geom_jitter(height = 0, size = 0.5, width = 0.05, aes(colour = scaledTime)) + 
  scale_colour_gradient(low = "#990000", high = "#FFCCCC") + ylim(42, 46)
p2

p1+p2


#Extended Data Fig. 3
data3 = read.csv(file = 'GC_subtype_avianFamily_data/avianIAV_persistent_GC_family.csv')
data3$family  = factor(data3$family, levels=c('Recurvirostridae',
                                              'Ciconiidae',
                                              'Threskiornithidae',
                                              'Phalacrocoracidae',
                                              'Scolopacidae',
                                              'Passeridae',
                                              'Alcidae',
                                              'Anatidae',
                                              'Struthionidae',
                                              'Numididae',
                                              'Psittaculidae',
                                              'Corvidae',
                                              'Accipitridae',
                                              'Sturnidae',
                                              'Ardeidae',
                                              'Podicipedidae',
                                              'Falconidae',
                                              'Columbidae',
                                              'Laridae',
                                              'Phasianidae',
                                              'Rallidae',
                                              'Gruidae',
                                              'ZZPersistent'))
p3 <- ggplot(data3, aes(factor(family), GC)) +
  geom_violin(width=1.5, size=0.4) + 
  geom_boxplot(width=0.1, size=0.5, color="grey", alpha=0.2, outliers = FALSE) +
  #geom_jitter(height = 0, size = 0.5,  width = 0.05, aes(colour = factor(Preidct))) +
  theme_classic() + ylim(42, 46) + theme(axis.text.x = element_text(angle = 90))
p3
