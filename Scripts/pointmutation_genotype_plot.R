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


## Part 1. Preparation for pastml, divide each seq to H5HP or nonH5HP based on defined HP in H5 large tree.
seglis=c('PB2','PB1','PA','NP','MP','NS')
for (seg in seglis ){
  namedf_loc<-read_tsv('H5HPclade230909.txt',col_names = T)
  genetaxa<-read_tsv(paste0('trait/',seg,'name.txt'),col_names=F)
  genetaxa<-merge(genetaxa,namedf_loc[,c("name","clade")],by.x="X1",by.y="name",all.x=T)
  genetaxa$subtype<-str_sub(str_split_fixed(as.character(genetaxa$X1),"\\|",5)[,3],5,)
  genetaxa$trait<-genetaxa$clade
  genetaxa[which(is.na(genetaxa$clade)),"trait"]='NonH5HP'
  write_tsv(genetaxa[,c("X1","trait")],paste0('trait/',seg,'_name_trait.txt'))
}

##run pastml 
#pastml
#pastml --tree nwktree/PB2_230909all.nwk --data trait/PB2_name_trait.txt --columns trait --prediction_method DOWNPASS
##script to get the transition nodes, define genotype group in each segment
#pastmltree.py 




## Part 2. process pastml results, and assign genotypes by grouping 8 segs
#combine segment number create genotype table
namedf_loc<-read_tsv('H5HPclade230909.txt',col_names = T)
sample<-filter(namedf_loc,grepl('^EPI',namedf_loc$name))
sample$name<-str_replace_all(sample$name,"\\'",'_')
sample$NAtype<-str_sub(str_split_fixed(sample$name,'\\|',5)[,3],3,4)
colnames(sample)<-c("name","HAclade","NAcladename")
sample$NAcladenum<-sample$NAcladename

for (seg in c('PB2','PB1','NP','MP','PA','NS')){
  seg_clade<-read_tsv(paste(seg,"_230909all_pastml/named.tree_",seg,"_230909all.nwkcladev2.txt",sep=""))
  colnames(seg_clade)<-c("name","origcladename")
  seg_clade$cladename<-paste(str_split_fixed(seg_clade$origcladename,'//',4)[,1],str_split_fixed(seg_clade$origcladename,'//',4)[,2],str_split_fixed(seg_clade$origcladename,'//',4)[,4],sep='//')
  segcladenum<-as.data.frame(cbind(unique(seg_clade$cladename),seq(1,length(unique(seg_clade$cladename)),1)))
  colnames(segcladenum)<-c('cladename','cladenum')
  segcladenum$cladenum<-paste(seg,segcladenum$cladenum,sep='_')
  write_tsv(segcladenum,paste(seg,'cladenum.txt',sep=""))
  seg_clade<-merge(seg_clade[,c("name","cladename")],segcladenum,by='cladename',all.x=T)
  sample<-merge(sample,seg_clade[,c('name','cladename','cladenum')],by='name',all.x=T)
  colnames(sample)[ncol(sample)]=paste(seg,'cladenum',sep = "")
  colnames(sample)[ncol(sample)-1]=paste(seg,'cladename',sep="")
}

write_tsv(sample,'Sample2344b_clade8segment.txt')
sample_rmna=drop_na(sample)
write_tsv(sample_rmna,'Sample2344b_clade8segment_rmna.txt')
groupsum=aggregate(sample_rmna$name, by=list(sample_rmna$HAclade,sample_rmna$PB2cladename,sample_rmna$PB2cladenum,sample_rmna$PB1cladename,sample_rmna$PB1cladenum,sample_rmna$PAcladename,sample_rmna$PAcladenum,sample_rmna$NPcladename,sample_rmna$NPcladenum,sample_rmna$NAcladename,sample_rmna$NAcladenum,sample_rmna$MPcladename,sample_rmna$MPcladenum,sample_rmna$NScladename,sample_rmna$NScladenum),length)

colnames(groupsum)<-c('HAclade','PB2cladename','PB2cladenum','PB1cladename','PB1cladenum','PAcladename','PAcladenum','NPcladename','NPcladenum','NAcladename','NAcladenum','MPcladename','MPcladenum','NScladename','NScladenum','x')
groupsum$genotype<-paste('G',seq(1,nrow(groupsum),1),sep='')
write_tsv(groupsum,'Sample2344b_groupsumv2.txt')
genotype=merge(groupsum,sample_rmna,by=c('HAclade','PB2cladenum','PB1cladenum','PAcladenum','NPcladenum','NAcladenum','MPcladenum','NScladenum'),all.y=T)
write_tsv(genotype,'Sample2344b_genotype.txt')
write_tsv(genotype[,c("name","HAclade","genotype")],'Sample2344b_genotype_simple.txt')
#count number of seqs in each genotype
ggplot(aes(x=genotype,y=log(x),fill=NAcladenum),data=groupsum)+geom_bar(stat='identity')+
  scale_y_continuous(expand = c(0,0))+theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())





## Part 3. Plot genotype
##NOTE: the node number is only applicable to the specific tree file, using for add node bootstrap for some important nodes. Please change the number for your tree.
#color
P32= createPalette(43,  c("#FFFFCC", "#FFCC00", "#CC3300","#CCFF00","#99CCFF","#33FFFF","#CCCCFF","#9900FF"))
P_mypal<-createPalette(28,c("#FED9A6","#B3E2CD","#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E" ,"#E6AB02" ,"#A6761D", "#666666"))
mypal=c("#FED9A6","#B3E2CD","#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E" ,"#E6AB02" ,"#A6761D", "#666666")
mypal2<-c("#f4bb8f", "#2d7a2c", "#09c553", "#86e6ca", "#1d686e", "#9ecbf4", "#6e31ac", "#c17aac", "#fa79f5", "#863563", "#427ff5","#EFE7BC", "#fcd107","#f7c5f1", "#3f16f9", "#add51f", "#6d4c2b", "#fe8f06", "#b34a01", "#e6262f", "#cc156d", "#fa1bfc", "#454366")
top9_mypal2=c("#8DD3C7",mypal2[1:10], "#D9D9D9", mypal2[11:20],"#BEBADA",P32[1:10] ,"#FB8072",P32[11:12], "#80B1D3", "#FDB462", "#B3DE69" ,"#FCCDE5","#B15928")
top14=c("#FDB462","#1B9E77", "#D95F02", "#7570B3","#66A61E" , "#E7298A", "#D9D9D9","#BEBADA","#FB8072", "#80B1D3","#8DD3C7" , "#B3DE69" ,"#FCCDE5","#B15928")
top9<-c("#FDB462", "#D9D9D9", "#BEBADA", "#FB8072", "#80B1D3", "#8DD3C7", "#B3DE69" ,"#FCCDE5","#B15928")
segcolor14<-c(brewer.pal(12,'Paired'),"#BABABA","#4D4D4D")



## Part 3.1.  large genotype heatmap plot for 1.5k seq tree
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

H5bs<-treeio::read.beast("sorted_H5_smallref_recent0909_1.5k.nexus")
trebsdf<-as_tibble(H5bs)

plabel<-ggtree(H5bs)%<+%GTdataname+geom_tiplab(aes(label=label_name),size=1)+geom_treescale(x=0,y=500,width=0.005)+
  geom_nodelab(aes(subset=(node%in%c(1521,1519,1518,1517,1516,1909,1955,1954,2012,2030,2049,1953,1952,1951,1515,1950)),x=branch,label=str_split_fixed(bb,'/',2)[,1]),size=1,vjust=-0.8)

p2<-gheatmap(plabel,genotype_test,offset=.01,width=.8,colnames=TRUE)+theme(legend.position="none")+
  scale_fill_manual(values=c(top9,unname(mypal2),"#FFFFCC", "#FFCC00", "#CC3300","#CCFF00","#99CCFF","#33FFFF","#CCCCFF","#9900FF"),na.value = 'white')
p
#save
ggsave('Genotype.pdf',width=6,height=55,limitsize = FALSE)



## Part3.2. 2k seq tree and host/prediction/country or other factors
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
  geom_tiplab(aes(subset=grepl('\\{',label_name),color=group,label=''),align = T,linetype = 1)+scale_color_manual(values=top14)+
  geom_tippoint(aes(subset=(label %in% pos$X1)),size=.5,color='red')+geom_treescale(x=0,y=500)+
  geom_nodelab(aes(subset=(node%in%c(2714,2715,3056,2711,3335,3363,2709,2708,3383,3382,2707,3823,3385,3386,3584,3656,3707,3773,2703,2701,2700,2698,2697,2695,2689,2687,2686)),x=branch,label=str_split_fixed(bb,'/',2)[,1]),size=1,vjust=-0.8)+geom_nodepoint(aes(subset=(node==2715)),color='red',shape=1)
p1label

hostp%>%insert_left(p1label,width=40)%>%insert_right(predp,width=1)

ggsave('Largetree.pdf',width=20,height=30)


##Part 3.3  Large tree annotating alrt value for each segment
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


#Part 3.4 plot mutation site
#Check mutation site comparing to reference
read_fas2df<-function(filename){
  df<-read_tsv(filename,col_names=F)
  odd=seq(1,nrow(df),2)
  even=seq(2,nrow(df),2)
  df<-data.frame(df[odd,1],df[even,1])
  print(table(str_sub(df$X1,1,1)))
  df$X1=str_sub(df$X1,2,)
  return(df)
}

GTgeno<-read_csv('GTtypefinal.csv',col_types = cols(.default = "c"))
sample<-read_tsv('Sample2344b_genotype0913v2_newGT.txt')
sample<-sample[,'name']
sample$EPI<-str_split_fixed(sample$name,'\\|',2)[,1]
GTgeno<-merge(GTgeno,sample,by.x='isolate_id',by.y='EPI',all.x=T)
GTgeno<-GTgeno%>%rename(label=name,Newgenotype=genotypenew)
mutlist<-read_tsv('Mutationlist.txt')

#read amino acid alignment and compare to Mutationlist
fas<-read_fas2df('sequence/H5_nt_230909align.translated.fas')
fas$EPI<-str_split_fixed(fas$X1,'\\|',2)[,1]
fas<-fas%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('HA',mutlist$col_name)]){
  fas[,i]<-str_sub(fas$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fas[,i])),table(fas[,i]),sep=':'),','),' '))
  }


faspb2<-read_fas2df('sequence/PB2_230909all.translated.fas')
faspb2$EPI<-str_split_fixed(faspb2$X1,'\\|',2)[,1]
faspb2<-faspb2%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PB2',mutlist$col_name)]){
  faspb2[,i]<-str_sub(faspb2$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspb2[,i])),table(faspb2[,i]),sep=':'),','),' '))
  }



faspb1<-read_fas2df('sequence/PB1_230909all.translated.fas')
faspb1$EPI<-str_split_fixed(faspb1$X1,'\\|',2)[,1]
faspb1<-faspb1%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PB1_',mutlist$col_name)]){
  faspb1[,i]<-str_sub(faspb1$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspb1[,i])),table(faspb1[,i]),sep=':'),','),' '))
}

faspb1f2<-read_fas2df('sequence/PB1-F2_230909align.translate.fas')
faspb1f2$EPI<-str_split_fixed(faspb1f2$X1,'\\|',2)[,1]
faspb1f2<-faspb1f2%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PB1-F2',mutlist$col_name)]){
  faspb1f2[,i]<-str_sub(faspb1f2$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspb1f2[,i])),table(faspb1f2[,i]),sep=':'),','),' '))
}


faspa<-read_fas2df('sequence/PA_230909all.translated.fas')
faspa$EPI<-str_split_fixed(faspa$X1,'\\|',2)[,1]
faspa<-faspa%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PA_',mutlist$col_name)]){
  faspa[,i]<-str_sub(faspa$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspa[,i])),table(faspa[,i]),sep=':'),','),' '))
}

fasnp<-read_fas2df('sequence/NP_230909all.translated.fas')
fasnp$EPI<-str_split_fixed(fasnp$X1,'\\|',2)[,1]
fasnp<-fasnp%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('NP',mutlist$col_name)]){
  fasnp[,i]<-str_sub(fasnp$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasnp[,i])),table(fasnp[,i]),sep=':'),','),' '))
}

fasna<-read_fas2df('sequence/N1_comb_230207_align0909edit.translated.fas')
fasna$EPI<-str_split_fixed(fasna$X1,'\\|',2)[,1]
fasna<-fasna%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('NA',mutlist$col_name)]){
  fasna[,i]<-str_sub(fasna$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasna[,i])),table(fasna[,i]),sep=':'),','),' '))
}

fasmp<-read_fas2df('sequence/MP_230909all.translated.fas')
fasmp$EPI<-str_split_fixed(fasmp$X1,'\\|',2)[,1]
fasmp<-fasmp%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('M1',mutlist$col_name)]){
  fasmp[,i]<-str_sub(fasmp$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasmp[,i])),table(fasmp[,i]),sep=':'),','),' '))
}

fasns<-read_fas2df('sequence/NS_230909all.translated.fas')
fasns$EPI<-str_split_fixed(fasns$X1,'\\|',2)[,1]
fasns<-fasns%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('NS1',mutlist$col_name)]){
  fasns[,i]<-str_sub(fasns$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasns[,i])),table(fasns[,i]),sep=':'),','),' '))
}

faspax<-read_fas2df('sequence/PA-X_230909all.translated.fas')
faspax$EPI<-str_split_fixed(faspax$X1,'\\|',2)[,1]
faspax<-faspax%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PA-X',mutlist$col_name)]){
  faspax[,i]<-str_sub(faspax$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspax[,i])),table(faspax[,i]),sep=':'),','),' '))
}

fasm2<-read_fas2df('sequence/M2_230909sample.translated.fas')
fasm2$EPI<-str_split_fixed(fasm2$X1,'\\|',2)[,1]
fasm2<-fasm2%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('M2',mutlist$col_name)]){
  fasm2[,i]<-str_sub(fasm2$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasm2[,i])),table(fasm2[,i]),sep=':'),','),' '))
}

fasns2<-read_fas2df('sequence/NS2_230909sample_align.translated.fas')
fasns2$EPI<-str_split_fixed(fasns2$X1,'\\|',2)[,1]
fasns2<-fasns2%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('NS2',mutlist$col_name)]){
  fasns2[,i]<-str_sub(fasns2$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasns2[,i])),table(fasns2[,i]),sep=':'),','),' '))
}


#write table
dfselect<-merge(sample,faspb2%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,faspb1%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,faspb1f2%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,faspa%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,faspax%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,fas%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,fasnp%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,fasna%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,fasmp%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,fasm2%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,fasns%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,fasns2%>%select(-c(X1,X1.1)),by='EPI',all.x=T)
dfselect<-merge(dfselect,GTgeno[,c('label','Newgenotype')],by.x='name',by.y='label',all.x=T)
write_xlsx(dfselect,'MutationBySeq1020.xlsx')


#plot mutation freq plot for each genotype
longdfselect<-melt(dfselect,id.vars = c('Newgenotype','name','EPI'),variable.name='variable',values.name='value')
longdfselect<-drop_na(longdfselect)
dfselect_mut_sum<-longdfselect%>%group_by(Newgenotype,variable,value)%>%summarise((n=n()))
mutlist<-read_tsv('Mutationlist.txt')
mutlist$order<-nrow(mutlist):1
mutlist$order<-as.numeric(mutlist$order)
mutlist$mutname<-paste0(mutlist$col_name,str_sub(unlist(str_match(mutlist$func,'[A-Z]\\d+?[A-Z]')),-1,-1))
for (k in which(grepl('/',unlist(str_split_fixed(mutlist$func,' ',2)[,1])))){
  print(k)
  m=str_count(unlist(str_match_all(mutlist$func[k],'\\d.+/[A-Z]')),'/')
  for (i in 1:m){
    add=mutlist[k,]
    add$mutname<-paste0(add$col_name,str_sub(unlist(str_split_fixed(add$func,' ',2)[,1]),-(2*i-1),-(2*i-1)))
    mutlist<-rbind(mutlist,add)
  }
}

mutlist$mutlabel<-paste(str_split_fixed(mutlist$col_name,'_',2)[,1],str_split_fixed(mutlist$func,' ',2)[,1],sep='_')
dfselect_mut_sum$mutname<-paste0(dfselect_mut_sum$variable,dfselect_mut_sum$value)
dfselect_select_mut<-filter(dfselect_mut_sum,mutname %in% c(mutlist$mutname))
colnames(dfselect_select_mut)[4]='num'
GTnum<-as.data.frame(table(dfselect$Newgenotype,exclude=NA))
colnames(GTnum)=c('Newgenotype','total')
dfselect_select_mut<-merge(dfselect_select_mut,GTnum,by='Newgenotype',all.x=T)
dfselect_select_mut<-merge(dfselect_select_mut,mutlist[,c('mutname','mutlabel','order')],by='mutname',all=T)
dfselect_select_mut$num[which(is.na(dfselect_select_mut$num))]<-0
dfselect_select_mut$Newgenotype[which(is.na(dfselect_select_mut$Newgenotype))]<-'G14'
dfselect_select_mut$total[which(is.na(dfselect_select_mut$total))]<-1
dfsum<-aggregate(num~Newgenotype+mutlabel+total+order,dfselect_select_mut,sum)
dfsum$prop<-dfsum$num/dfsum$total
dfsum$Newgenotype<-factor(dfsum$Newgenotype,levels=c(paste('G',1:14,sep='')))
ggplot(dfsum,aes(x=Newgenotype,y=reorder(mutlabel,order),fill=prop))+
  geom_tile(width=0.9)+theme_bw()+
  scale_fill_gradient(low="white",high="#4D4D4D",
                      guide = guide_colorbar(label = TRUE,frame.colour = "black",ticks = TRUE,ticks.colour = 'black'))+
  theme(panel.grid = element_blank())+scale_x_discrete(expand = c(0, 0))+xlab('')+ylab('')

ggsave('mutation.pdf',height=8,width=6)
