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
fas<-read_fas2df('sequence\\H5_nt_230909align.translated.fas')
fas$EPI<-str_split_fixed(fas$X1,'\\|',2)[,1]
fas<-fas%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('HA',mutlist$col_name)]){
  fas[,i]<-str_sub(fas$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fas[,i])),table(fas[,i]),sep=':'),','),' '))
  }


faspb2<-read_fas2df('sequence\\PB2_230909all.translated.fas')
faspb2$EPI<-str_split_fixed(faspb2$X1,'\\|',2)[,1]
faspb2<-faspb2%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PB2',mutlist$col_name)]){
  faspb2[,i]<-str_sub(faspb2$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspb2[,i])),table(faspb2[,i]),sep=':'),','),' '))
  }



faspb1<-read_fas2df('sequence\\PB1_230909all.translated.fas')
faspb1$EPI<-str_split_fixed(faspb1$X1,'\\|',2)[,1]
faspb1<-faspb1%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PB1_',mutlist$col_name)]){
  faspb1[,i]<-str_sub(faspb1$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspb1[,i])),table(faspb1[,i]),sep=':'),','),' '))
}

faspb1f2<-read_fas2df('sequence\\PB1-F2_230909align.translate.fas')
faspb1f2$EPI<-str_split_fixed(faspb1f2$X1,'\\|',2)[,1]
faspb1f2<-faspb1f2%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PB1-F2',mutlist$col_name)]){
  faspb1f2[,i]<-str_sub(faspb1f2$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspb1f2[,i])),table(faspb1f2[,i]),sep=':'),','),' '))
}


faspa<-read_fas2df('sequence\\PA_230909all.translated.fas')
faspa$EPI<-str_split_fixed(faspa$X1,'\\|',2)[,1]
faspa<-faspa%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('PA_',mutlist$col_name)]){
  faspa[,i]<-str_sub(faspa$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(faspa[,i])),table(faspa[,i]),sep=':'),','),' '))
}

fasnp<-read_fas2df('sequence\\NP_230909all.translated.fas')
fasnp$EPI<-str_split_fixed(fasnp$X1,'\\|',2)[,1]
fasnp<-fasnp%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('NP',mutlist$col_name)]){
  fasnp[,i]<-str_sub(fasnp$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasnp[,i])),table(fasnp[,i]),sep=':'),','),' '))
}

fasna<-read_fas2df('sequence\\N1_comb_230207_align0909edit.translated.fas')
fasna$EPI<-str_split_fixed(fasna$X1,'\\|',2)[,1]
fasna<-fasna%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('NA',mutlist$col_name)]){
  fasna[,i]<-str_sub(fasna$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasna[,i])),table(fasna[,i]),sep=':'),','),' '))
}

fasmp<-read_fas2df('sequence\\MP_230909all.translated.fas')
fasmp$EPI<-str_split_fixed(fasmp$X1,'\\|',2)[,1]
fasmp<-fasmp%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('M1',mutlist$col_name)]){
  fasmp[,i]<-str_sub(fasmp$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasmp[,i])),table(fasmp[,i]),sep=':'),','),' '))
}

fasns<-read_fas2df('sequence\\NS_230909all.translated.fas')
fasns$EPI<-str_split_fixed(fasns$X1,'\\|',2)[,1]
fasns<-fasns%>%dplyr::filter(EPI %in% unlist(sample[,'EPI']))
for (i in mutlist$col_name[grepl('NS1',mutlist$col_name)]){
  fasns[,i]<-str_sub(fasns$X1.1,mutlist$site[which(mutlist$col_name==i)],mutlist$site[which(mutlist$col_name==i)])
  print(paste(i,str_flatten(paste(names(table(fasns[,i])),table(fasns[,i]),sep=':'),','),' '))
}

faspax<-read_fas2df('sequence\\PA-X_230909all.translated.fas')
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

fasns2<-read_fas2df('sequence\\NS2_230909sample_align.translated.fas')
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
