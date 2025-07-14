library(sf)
library(sp)
library(maptools)
library(ggmap)
library(cowplot)

#read geojson format map
xtemp1<-st_read('Europelevel2.geo.json')
xtemp2<-st_read('Africa.geo.json')
xtemp3<-st_read('Asia.geo.json')
xtemp<-rbind(xtemp1,xtemp2,xtemp3)
ggplot(xtemp)+geom_sf(aes(geometry=geometry),size=0.1)+coord_sf(xlim=c(-40,55),ylim=c(30,80))+scale_x_continuous(expand=c(0,0))

#read bayes factor calculated from beast results
bf<-read_tsv('beast/bf123.txt')
colnames(bf)=c('from','to','bayes-factor','postprob')
rate<-read_tsv('beast/H5G2rate.txt')
coord<-read_tsv('coord_df.tsv',col_names = FALSE)
colnames(coord)=c('location','lat','lon','adm0_a3')
countrynum<-read_tsv('countrynum.txt')
ratedf<-merge(rate,bf,by=c('from','to'),all.x=T)
ratedf<-merge(ratedf,coord,by.x='from',by.y='location',all.x=T)
ratedf<-merge(ratedf,coord,by.x='to',by.y='location',all.x=T)
ratedf<-ratedf%>%rename(bf=`bayes-factor`)
#select BF >10 or >30
ratedf10<-filter(ratedf,bf>10)
ratedf10$bfcut<-cut(ratedf10$bf,breaks=c(10,30,100,1000,1000000),labels = c('10-30','30-100','100-1000','>1000'))
ratedf30<-ratedf10%>%filter(bf>30)

#plot map
theme_set(theme_bw())
xtempmap<-merge(xtemp,countrynum,by='adm0_a3',all.x=T)
ggplot(xtempmap)+
  geom_sf(aes(geometry=geometry,fill=num),size=0.1)+
  coord_sf(xlim=c(-20,45),ylim=c(10,70))+
  scale_x_continuous(expand=c(0,0))+
  geom_curve(data=ratedf10,
             aes(x=lon.x,y=lat.x,xend=lon.y,yend=lat.y,color=bfcut,size=median),
             arrow=arrow(length=unit(0.15,'cm')),
             curvature=0.2)+
  scale_color_manual(values=
                       c( "#6BAED6" ,"#4292C6", "#2171B5","#08519C","#08306B"),name='Bayes factor'
  )+
  scale_size(range=c(0.1,0.8),name='Median transmission rate')+
  xlab('')+ylab('')+
  scale_fill_continuous(low='#FFFFCC',high="#FC4E2A",guide='colorbar',na.value='#F0F0F0',name='Number of sequence')
#save
ggsave('beast/mapbf0920.pdf',width = 6,height = 10)