library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(lubridate)
library(ggtree)
library(ggplot2)
library(treeio)
library(reshape2)
library(Polychrome)
library(RColorBrewer)
library(aplot)

#color
top14=c("#FDB462","#1B9E77", "#D95F02", "#7570B3","#66A61E" , "#E7298A", "#D9D9D9","#BEBADA","#FB8072", "#80B1D3","#8DD3C7" , "#B3DE69" ,"#FCCDE5","#B15928")
segcolor14<-c(brewer.pal(12,'Paired'),"#BABABA","#4D4D4D")

##Part 1: Extended Fig. 6b
## large genotype heatmap plot for 1.5k seq tree
# This file 'Sample2344b_genotype0913v2_newGT.txt' is genotype by my definition, only use below for select target samples. 
# GTtypefinal.csv is the final genotype
sample<-read_tsv('Sample2344b_genotype0913v2_newGT.txt')
sample<-sample[,'name']
GTdata<-read_csv('GTtypefinal.csv',col_types = cols(.default = "c"))
sample$EPI<-str_split_fixed(sample$name,'\\|',2)[,1]
GTdata<-merge(GTdata,sample,by.x='isolate_id',by.y='EPI',all.x=T)
GTdata<-GTdata%>%rename(label=name,Newgenotype=genotypenew)
GTdata$label_name<-paste(str_split_fixed(GTdata$label,'\\|',6)[,6],'_{',GTdata$Newgenotype,'}',sep='')
GTdata$label_name[which(is.na(GTdata$Newgenotype))]=str_split_fixed(GTdata$label[which(is.na(GTdata$Newgenotype))],'\\|',6)[,6]
GTdataname<-GTdata[,c("label","label_name")]
genotype_test<-GTdata%>%select(label,PB2:NS)
rownames(genotype_test)=genotype_test[,1]
genotype_test=genotype_test[,-1]

##  2k seq tree and host/prediction/country or other factors
#host prediction heatmap
namedf_loc<-read_tsv('H5HPclade230909.txt',col_names = T)
metadata1<-readxl::read_xlsx('seqname_host_prediction.xlsx',sheet=1)
metadata2<-readxl::read_xlsx('seqname_host_prediction.xlsx',sheet=2)
metadata<-rbind(metadata1,metadata2)
metadata$seqname<-str_replace_all(metadata$seqname,'>','')
metadata$Prediction<-as.character(metadata$Prediction)
metadatafig<-merge(metadata,namedf_loc,by.x='seqname',by.y='name',all.x=T)
metadatafig<-filter(metadatafig,clade=='2.3.4.4b')
metadatafig$Prediction[which(metadatafig$Prediction==-1)]<-'Negative'
metadatafig$Prediction[which(metadatafig$Prediction==1)]<-'Positive'

H5bs2<-treeio::read.beast("sorted_H5_smallref_recent0909_0906.nexus")
trebsdf2<-as_tibble(H5bs2)

seqname<-as.data.frame(H5bs2@phylo$tip.label)
colnames(seqname)=c('seqname')
metadatafig<-merge(metadatafig,seqname,by.x='seqname',by.y='seqname',all=T)
colnames(metadatafig)[1]='label'
hostp<-ggplot(aes(x=1,y=label,fill=Host),data=metadatafig)+
  geom_tile()+ylab('')+xlab('Host')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+ 
  scale_fill_discrete(breaks=c('Fox', 'Gull', 'Mink', 'Human', 'OtherMammal','OtherAvian'))+
  scale_fill_manual(values=c("Fox"="#FF0000","Gull"="#A9A9A9","Mink"="#FF3333","Human"="#CC0000",'OtherMammal'="#FF6666",'OtherAvian'="#C0C0C0"),na.value='white')

pos<-read_table('firstpositive.txt',col_names = FALSE)
predp<-ggplot(aes(x=1,y=label,fill=Prediction),data=metadatafig)+
  geom_tile()+ylab('')+xlab('Prediction')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_discrete(breaks=c('Positive','Negative'))+scale_fill_manual(values=c("Positive"="#E7211A" ,"Negative"="#D0D1D1"),na.value='white')

groupInfo<- split(GTdata$label, GTdata$Newgenotype)
H5bs2<-groupOTU(H5bs2,groupInfo)
GTdata<-GTdata%>%select(label,everything())
p1label<-ggtree(H5bs2)%<+%GTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1)+scale_color_manual(values=top14)+geom_tippoint(aes(subset=(label %in% pos$X1)),size=.5,color='red')+geom_treescale(x=0,y=500) + 
  geom_nodelab(aes(subset=(node%in%c(2714,2715,3056,2711,3335,3363,2709,2708,3383,3382,2707,3823,3385,3386,3584,3656,3707,3773,2703,2701,2700,2698,2697,2695,2689,2687,2686)),x=branch,label=str_split_fixed(bb,'/',2)[,1]),size=1,vjust=-0.8)+geom_nodepoint(aes(subset=(node==2715)),color='red',shape=1)

locationp<-ggplot(aes(x=1,y=label,fill=Location),data=metadatafig)+
  geom_tile()+ylab('')+xlab('Location')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+ 
  scale_fill_discrete(breaks=c('Europe', 'Africa', 'Asia', 'NorthAmerica', 'SouthAmerica'))+
  scale_fill_manual(values=c("Europe"="#D95F02","Africa"="#FB8072","Asia"="#FDB462","NorthAmerica"="#80B1D3",'SouthAmerica'="#BEBADA"),na.value='white')


hostp%>%insert_left(p1label,width=40)%>%insert_right(predp,width=1)%>%insert_right(locationp,width=2)

ggsave('Largetree.pdf',width=20,height=30)


##Part 2  Large tree annotating alrt value for each segment, Supplementary Fig. 16-23
#prepare Genotype group to color the tree tips by genotype
sample<-read_tsv('Sample2344b_genotype0913v2_newGT.txt')
sample<-sample[,'name']
GTdata<-read_csv('GTtypefinal.csv',col_types = cols(.default = "c"))
sample$EPI<-str_split_fixed(sample$name,'\\|',2)[,1]
GTdata<-merge(GTdata,sample,by.x='isolate_id',by.y='EPI',all.x=T)
segGTdata<-GTdata%>%rename(label=name,Newgenotype=genotypenew)
segGTdata$label_name<-paste(str_split_fixed(segGTdata$label,'\\|',6)[,6],'_{',segGTdata$Newgenotype,'}',sep='')
segGTdata<-segGTdata%>%select(label,everything())
segGT<-segGTdata[,c('label','Newgenotype')]
segGT<-segGT%>%drop_na(Newgenotype)
segGT$label_name<-paste(str_split_fixed(segGT$label,'\\|',6)[,6],'_{',segGT$Newgenotype,'}',sep='')
groupInfoseg <- split(segGT$label, segGT$Newgenotype)

#segment tree v2 PB2 
segtree<-treeio::read.beast('/alrtTree/PB2_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','PB2')],by='label',all.x=T)
segbardata$PB2[which(is.na(segbardata$Newgenotype))]=NA
segbar<-ggplot(aes(x=1,y=label,fill=PB2),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(41553,46309,41520,41519,41516,41502,66542,66551,66539)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('PB2segv2.pdf',width=15,height=50,limitsize = FALSE)


#segment tree v2 PB1
segtree<-treeio::read.beast('/alrtTree/PB1_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','PB1')],by='label',all.x=T)
segbardata$PB1[which(is.na(segbardata$Newgenotype))]=NA
segbar<-ggplot(aes(x=1,y=label,fill=PB1),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

ggtree(segtree)+geom_nodelab(aes(x=branch,label=node),size=2,hjust=.5)
#ggsave('PB1nodelab.pdf',width=15,height=50,limitsize = FALSE) #view node label

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(41423,41424,41424,41425,41426,41427,41428,41429,41433,41434,41422,41421,58606,58605,66530,75987,66537,66604)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('PB1segv2.pdf',width=15,height=50,limitsize = FALSE)


#segment tree v2 PA
segtree<-treeio::read.beast('/alrtTree/PA_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','PA')],by='label',all.x=T)
segbardata$PA[which(is.na(segbardata$Newgenotype))]=NA
segbar<-ggplot(aes(x=1,y=label,fill=PA),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

ggtree(segtree)+geom_nodelab(aes(x=branch,label=node),size=2,hjust=.5)
#ggsave('PAnodelab.pdf',width=15,height=50,limitsize = FALSE)

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(41625,41626,41627,43930,47925,41576,41577,41580,41578,76704,76729,59924)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('PAsegv2.pdf',width=15,height=50,limitsize = FALSE)


#segment tree v2 NP
segtree<-treeio::read.beast('/alrtTree/NP_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','NP')],by='label',all.x=T)
segbardata$NP[which(is.na(segbardata$Newgenotype))]=NA
segbar<-ggplot(aes(x=1,y=label,fill=NP),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

ggtree(segtree)+geom_nodelab(aes(x=branch,label=node),size=2,hjust=.5)
#ggsave('NPnodelab.pdf',width=15,height=50,limitsize = FALSE)

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(42008,42009,42010,42011,42012,42013,42016,42017,65935,75513,75514,65997,65996,65941,65938,53734,55481,53735,53736)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('NPsegv2.pdf',width=15,height=50,limitsize = FALSE)


#segment tree v2 MP
segtree<-treeio::read.beast('/alrtTree/MP_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','MP')],by='label',all.x=T)
segbardata$MP[which(is.na(segbardata$Newgenotype))]=NA
segbar<-ggplot(aes(x=1,y=label,fill=MP),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

ggtree(segtree)+geom_nodelab(aes(x=branch,label=node),size=2,hjust=.5)
#ggsave('MPnodelab.pdf',width=15,height=50,limitsize = FALSE)

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(43147,43148,43149,43153,67747,67723,43154,58113,58118,43162,79343,79389,79390)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('MPsegv2.pdf',width=15,height=50,limitsize = FALSE)



#segment tree v2 NS
segtree<-treeio::read.beast('/alrtTree/NS_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','NS')],by='label',all.x=T)
segbardata$NS[which(is.na(segbardata$Newgenotype))]=NA
segbar<-ggplot(aes(x=1,y=label,fill=NS),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

ggtree(segtree)+geom_nodelab(aes(x=branch,label=node),size=1,hjust=1)
#ggsave('NSnodelab.pdf',width=15,height=50,limitsize = FALSE)

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(42771,42797,62913,62911,49654,64277,64276,49651,84711,84172,84173)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('NSsegv2.pdf',width=15,height=50,limitsize = FALSE)



#segment tree v2 N1
segtree<-treeio::read.beast('/alrtTree/N1_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','NA')],by='label',all.x=T)

segbardata[,'NA'][which(is.na(segbardata$Newgenotype))]=NA
colnames(segbardata)[3]='N1'
segbar<-ggplot(aes(x=1,y=label,fill=N1),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

ggtree(segtree)+geom_nodelab(aes(x=branch,label=node),size=1,hjust=1)
#ggsave('N1nodelab.pdf',width=15,height=50,limitsize = FALSE)

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(15556,15559,15560,21080,18622,15567,15585,15566,15565)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('N1segv2.pdf',width=15,height=50,limitsize = FALSE)



#segment tree v2 H5
segtree<-treeio::read.beast('/alrtTree/H5_230909all.nexus')
segtreedf<-as_tibble(segtree)
segtree<-groupOTU(segtree,groupInfoseg)

segbardata<-merge(data.frame(label=segtree@phylo$tip.label),segGTdata[,c('label','Newgenotype','HA')],by='label',all.x=T)

segbardata[,'HA'][which(is.na(segbardata$Newgenotype))]=NA
segbar<-ggplot(aes(x=1,y=label,fill=HA),data=segbardata)+
  geom_tile()+ylab('')+xlab('')+
  theme(axis.text= element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=unname(segcolor14),na.value='white')

ggtree(segtree)+geom_nodelab(aes(x=branch,label=node),size=1,hjust=1)
#ggsave('H5nodelab.pdf',width=15,height=50,limitsize = FALSE)

pGT<-ggtree(segtree)%<+%segGTdata+
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1,linesize = .001)+
  scale_color_manual(values=c('white',top14),na.value='white')+
  geom_nodelab(aes(subset=(node%in%c(21379,21385,32738,21404,26161,21425,21450,27270,21453)),x=branch,label=bb),size=2,vjust=-0.8)+
  geom_treescale(x=0,y=500,width=0.05)#+theme(legend.position = 'none')
pGT

segbar%>%insert_left(pGT,width=10)
ggsave('H5segv2.pdf',width=15,height=50,limitsize = FALSE)