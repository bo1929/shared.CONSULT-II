require(reshape2); require(ggplot2)


dtc=read.csv('dist-to-closest.txt',sep="\t",he=F)
dtc$V1 = sub(".fna","",dtc$V1)
tax=read.csv('rank_queries.tsv',sep="\t")
dtc = merge(dtc,tax,by.x="V1",by.y="genome")
nrow(dtc)

summary(cut(dtc$V3,c(0,0.001,0.02,0.06,0.12,0.16,0.22,0.35),right=F))
summary(cut(dtc$V3,c(-0.000001,0.0001,0.02,0.06,0.12,0.16,0.2,0.25,0.5)))

qplot(reorder(sub(".fna","",V1),V3),V3,data=dtc)+theme_classic()+
  geom_hline(yintercept=c(0.0001,0.02,0.06,0.12,0.16,0.2,0.25),color="red",linetype=2)+
  theme(axis.text.x=element_text(angle=90))+xlab("qyery")+ylab("Mash distance to cloests")
ggsave("query-distances.pdf",width=12,height=6)



head(tax)
qplot(reorder(phylum,phylum,
              function(x)-length(x)),fill=as.factor(iorder ),data=tax)+
  geom_bar(color="black")+
  theme_classic()+
  theme(legend.position = "none",axis.text.x=element_text(angle=90,hjust = 1))+
  xlab("")
ggsave("query-groups.pdf",width = 5,height = 6)

unique(tax$phylum)

head(dtc)
k=read.csv('kraken_results.csv')
k = dcast(data=melt(k[2:9],id.vars = 1:2),formula = variable+genome~value)
k$m= "Kraken-II"
head(k) 
c=read.csv('consult_results.csv')
c = dcast(data=melt(c[2:9],id.vars = 1:2),formula = variable+genome~value)
c$m= "CONSULT-II"
head(c)
b=rbind(k,c)
head(b)

k2 = merge(dtc[,1:13],b,by.x = "V1",by.y="genome")
head(k2)
names(k2)[14]="level"
names(k2)[5:13] = sub("^i","",names(k2)[5:13])
k2=melt(k2,measure.vars = c("FN","FP","TN","TP"))
head(k2)
k2$value = apply(k2,1,function(x) ifelse(as.numeric(x[[x[["level"]]]])==0,0,as.numeric(x[["value"]])))

#k = merge(k,tax[,1:10],by.x="V1",by.y="igenome")
#k2=k
#k2$phylum = ifelse(k2$iphylum == 0, NA, k$phylum)
#k2$class = ifelse(k$iclass == 0, NA, k$class)
#k2$order = ifelse(k$iorder == 0, NA, k$order)
#k2$family = ifelse(k$ifamily == 0, NA, k$family)
#k2$genus = ifelse(k$igenus == 0, NA, k$genus)
#head(k2)

nozeroos = dtc[apply(dtc[,8:13]==0,1,sum)==0,"V1"]

ks = dcast(data=k2[,c("level","m","variable","value","V3")],
           formula = cut(V3,c(0,0.001,0.02,0.06,0.12,0.16,0.22,0.35),right=F)+level+m~variable,
          value.var = "value",fun.aggregate = sum)
names(ks)[1]="bin"
ks$s = apply(ks,1,function(x) sum(as.numeric(x[4:7])))
ks

#ks=ks[ks$level!="kingdom",]
#ks$variable=factor(ks$variable,c("kingdom" ,"phylum", "class", "order" ,"family", "genus" ,"species"))


ggplot(aes(x=bin,y=TP/(TP+FN),
           color=m,linetype=m,shape=m),data=ks)+
  geom_point()+
  facet_wrap(~level)+
  geom_line(aes(group=interaction(m)))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_classic()+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme(legend.position = "bottom" ,
        axis.text.x = element_text(angle=90))+ # c(0.1,0.3))+
  xlab("Distance to closest")+ylab("Recall")
ggsave("recall.pdf",width=6.5,height = 5)


ggplot(aes(x=bin,y=TP/(TP+FP),
           color=m,linetype=m,shape=m),data=ks)+
  geom_point()+
  facet_wrap(~level)+
  geom_line(aes(group=interaction(m)))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_classic()+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme(legend.position = "bottom" ,
        axis.text.x = element_text(angle=90))+ # c(0.1,0.3))+
  xlab("Distance to closest")+ylab("Precision")
ggsave("precision.pdf",width=6.5,height = 5)

ggplot(aes(x=bin,y=2*TP/(2*TP+FP+FN),
           color=level,linetype=m,shape=m),data=ks)+
  geom_point()+
  facet_wrap(~floor((as.numeric(ks$level))/2.5),
            labeller = function(x) {x[1,]="Low";x[1,]="Medium";x[3,]="High"})+
  geom_line(aes(group=interaction(m,level)))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_classic()+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme(legend.position = "bottom" ,
        axis.text.x = element_text(angle=90,hjust=1),
        legend.box.margin = margin(0),
        legend.margin = margin(0),
        )+ # c(0.1,0.3))+
  xlab("Distance to closest")+ylab("F1")
ggsave("F1.pdf",width=7,height = 5)



ggplot(aes(y=TP/(TP+FN),x=TP/(TP+FP),color=bin,shape=level,linetype=m),data=ks)+
  stat_summary()+
  stat_summary(aes(group=interaction(bin,m)),geom="line")+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme_classic()+xlab("Precision")+ylab("Recall")+
  theme(legend.position = "bottom",legend.direction = "horizontal")
ggsave("precision-recall.pdf",width=7,height = 5)


ggplot(aes(y=TP/(TP+FN),x=TP/(TP+FP)),
       data=ks)+
  geom_line(aes(group=m,linetype=m),size=0.7)+
  geom_point(aes(color=bin,shape=m),size=2)+
  facet_wrap(~level)+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  scale_linetype(name="")+
  geom_line(aes(group=bin),color="grey50",alpha=0.5,size=0.2)+
  theme_classic()+xlab("Precision")+ylab("Recall")+
  theme(legend.position = "bottom",legend.direction = "horizontal")
ggsave("precision-recall-methods.pdf",width=7,height = 6)


ggplot(aes(y=TP/(TP+FN),x=TP/(TP+FP)), data=ks)+
  geom_line(aes(group=m,linetype=m),size=0.7)+
  geom_point(aes(color=level,shape=m),size=2)+
  facet_wrap(~bin)+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  scale_linetype(name="")+
  geom_line(aes(group=level),color="grey50",alpha=0.5,size=0.2)+
  theme_classic()+xlab("Precision")+ylab("Recall")+
  theme(legend.position = "bottom",legend.direction = "horizontal")
ggsave("precision-recall-methods2.pdf",width=7,height = 6)
