require(reshape2); require(ggplot2)


dtc=read.csv('dist-to-closest.txt',sep="\t",he=F)
dtc$V1 = sub(".fna","",dtc$V1)
tax=read.csv('rank_queries.tsv',sep="\t")
dtc = merge(dtc,tax,by.x="V1",by.y="genome")
nrow(dtc)

summary(cut(dtc$V3,c(0,0.001,0.02,0.06,0.12,0.16,0.22,0.35),right=F))
summary(cut(dtc$V3,c(-0.000001,0.0001,0.02,0.06,0.12,0.16,0.2,0.25,0.5)))

qplot(reorder(sub(".fna","",V1),V3),V3,data=dtc)+theme_classic()+
  geom_hline(yintercept=c(0.001,0.02,0.06,0.12,0.16,0.22),color="red",linetype=2)+
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
cl=read.csv('CLARK-bacteria-eval.csv')
cl = dcast(data=melt(cl[2:9],id.vars = 1:2),formula = variable+genome~value)
head(cl)
cl$m= "CLARK"
readcons = function(f,n) {
  c=read.csv(f)
  c = c[,c(1:3,6:12)]
  c = dcast(data=melt(c[2:9],id.vars = 1:2),formula = variable+genome~value)
  c$m = n
  c
}
c01 = readcons('CONSULT-bacteria-eval-th05_d1_c001.csv','CONSULT-II (0.01)')
c03 = readcons('CONSULT-bacteria-eval-th05_d1_c003.csv','CONSULT-II (0.03)')
c = readcons('CONSULT-bacteria-eval-th05_d1.csv','CONSULT-II (0.00)')
b=rbind(k,c01,c03,c,cl)
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
ks$m = factor(ks$m,levels=c("CONSULT-II (0.00)","CONSULT-II (0.01)","CONSULT-II (0.03)", "Kraken-II","CLARK"))
ks

#ks=ks[ks$level!="kingdom",]
#ks$variable=factor(ks$variable,c("kingdom" ,"phylum", "class", "order" ,"family", "genus" ,"species"))


ggplot(aes(x=bin,y=2*TP/(2*TP+FP+FN),
           color=level,linetype=m,shape=m),data=ks[grepl("CONS",ks$m),])+
  geom_point(alpha=0.7)+
  facet_wrap(~floor((as.numeric(level))/2.5),
        labeller = function(x) {x[1,1]="High resolution";x[2,1]="Medium resolution";x[3,1]="Low resolution";x})+
  geom_line(aes(group=interaction(m,level)))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_classic()+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme(legend.position = "bottom" ,
        axis.text.x = element_text(angle=90,hjust=1),
        legend.box.margin = margin(0),
        legend.margin = margin(0))+ 
  xlab("Distance to closest")+ylab("F1")
ggsave("F1-consults.pdf",width=7,height = 4.5)



ggplot(aes(x=bin,y=2*TP/(2*TP+FP+FN),
        color=sub(" .*","",m),shape=sub(" .*","",m),
        linetype=sub("^$","*",gsub("[^0-9.]",'',m))),
       data=ks)+
  facet_wrap(~level)+
  geom_point(alpha=0.7)+
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
ggsave("F1-all.pdf",width=7.3,height = 5.5)

ggplot(aes(y=TP/(TP+FN),x=TP/(TP+FP),
           color=sub(" .*","",m),#shape=sub(" .*","",m),
           shape=sub("^$","*",gsub("[^0-9.]",'',m))),
       data=ks)+
  facet_wrap(~level)+
  #geom_path(aes(group=m,linetype=m),size=0.3)+
  geom_path(aes(group=bin),color="grey50",linetype=1,alpha=0.5,size=0.3)+
  geom_point(aes(),size=2,alpha=0.75)+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme_classic()+xlab("Precision")+ylab("Recall")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.margin = margin(0),
        legend.margin = margin(0))+
  guides(linetype=guide_legend(nrow=1, byrow=TRUE),
         shape=guide_legend(nrow=1, byrow=TRUE))
ggsave("precision-recall-methods-all.pdf",width=7,height = 5)


ggplot(aes(y=TP/(TP+FN),x=TP/(TP+FP)),
       data=ks[grepl("CONS",ks$m),])+
  facet_wrap(~level)+
  #geom_path(aes(group=m,linetype=m),size=0.3)+
  geom_path(aes(group=bin),color="grey50",linetype=1,alpha=0.5,size=0.3)+
  geom_point(aes(color=bin,shape=gsub("CONSULT-II ","",m)),size=2,alpha=0.75)+
  scale_color_brewer(palette = "Paired",name="")+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme_classic()+xlab("Precision")+ylab("Recall")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.margin = margin(0),
        legend.margin = margin(0))+
  guides(linetype=guide_legend(nrow=2, byrow=TRUE),
         shape=guide_legend(nrow=2, byrow=TRUE))
ggsave("precision-recall-methods-consult.pdf",width=7,height = 5)


ggplot(aes(y=TP/(TP+FN),x=TP/(TP+FP)),
       data=ks[!grepl("0[01]",ks$m ),])+
  facet_wrap(~level)+
  #geom_path(aes(group=m,linetype=m),size=0.3)+
  geom_path(aes(group=bin),color="grey50",linetype=1,alpha=0.5,size=0.3)+
  geom_point(aes(color=bin,shape=m),size=2,alpha=0.75)+
  scale_color_brewer(palette = "Paired",name="")+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme_classic()+xlab("Precision")+ylab("Recall")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.margin = margin(0),
        legend.box.spacing = unit(0,"pt"),
        legend.margin = margin(0))+
  guides(linetype=guide_legend(nrow=2, byrow=TRUE),
         shape=guide_legend(nrow=2, byrow=TRUE))
ggsave("precision-recall-methods-03.pdf",width=7.3,height = 5)

ggplot(aes(y=TP/(TP+FN),x=TP/(TP+FP),color=bin,shape=level,linetype=m),data=ks)+
  stat_summary()+
  stat_summary(aes(group=interaction(bin,m)),geom="line")+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme_classic()+xlab("Precision")+ylab("Recall")+
  theme(legend.position = "bottom",legend.direction = "horizontal")
ggsave("precision-recall.pdf",width=7,height = 5)



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

ggplot(aes(x=V3,y=2*TP/(2*TP+FP+FN),
           color=level,
           linetype=factor(m,levels=c("CONSULT-II", "Kraken-II","CLARK")),
              shape=factor(m,levels=c("CONSULT-II", "Kraken-II","CLARK"))),
       data= dcast(data=k2[,c("level","m","variable","value","V3")], formula = V3+level+m~variable,
                   value.var = "value",fun.aggregate = sum))+
  geom_point(alpha=0.7)+
  facet_wrap(~level#, 
             #labeller = function(x) {x[1,1]="High resolution";x[2,1]="Medium resolution";x[3,1]="Low resolution";x}
    )+
  #geom_line(aes(group=interaction(m,level)))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_classic()+
  stat_smooth(se=F,span = 0.7,method="glm",
              method.args=list(family=binomial),size=0.75)+
  scale_shape(name="")+
  scale_linetype(name="")+
  theme(legend.position = "bottom" ,
        axis.text.x = element_text(angle=90,hjust=1),
        legend.box.margin = margin(0),
        legend.margin = margin(0),
  )+ # c(0.1,0.3))+
  xlab("Distance to closest")+ylab("F1")+
  guides(linetype=guide_legend(nrow=2, byrow=TRUE),
         shape=guide_legend(nrow=2, byrow=TRUE))
ggsave("F1-nobin.pdf",width=7,height = 5)

km = dcast(data=k2[,c("level","m","variable","value","V3")], formula = V3+level+m~variable,
      value.var = "value",fun.aggregate = sum)
head(km)
with(km, km[m=="Kraken-II"&V3<0.02&TP/(TP+FN) < 0.3 ,])


cvv = read.csv('CONSULT-bacteria-eval-th05_d1.csv',sep=",",he=T)
head(cvv)
cvv = melt(cvv,id.vars = 1:5)
cvv = cvv[cvv$variable!="superkingdom",]
cvvp = cvv[cvv$value%in% c("TP","FP"),]
cvvn = cvv[!cvv$value%in% c("TP","FP"),]
#dcast(data=cvv,formula=cut(V,5)+)
ggplot(aes(x=vote,color=variable,linetype=value),data=cvvp)+
  stat_ecdf()+
  #facet_wrap(~variable)+
  #geom_vline(xintercept = c(0.03,0.01),color="grey40",linetype=3)+
  theme_classic()+
  scale_y_continuous(name="ECDF")+
  scale_x_continuous(name="vote",trans="log10")+
  scale_linetype_manual(name="",values=c(1,3))+
  annotate("rect", xmin = 0.003, xmax = 0.01, ymin = 0, ymax = 1, alpha = .1)+
  annotate("rect", xmin = 0.003, xmax = 0.03, ymin = 0, ymax = 1, alpha = .1)+
  scale_color_brewer(palette = "Dark2",name="")
  #scale_color_manual(values=c("#AA3333","#22DD55","#DD6666","#225599"),name="")
ggsave("votes.pdf",width=7,height=5.5)

ggplot(aes(x=vote,color=variable,linetype=value),data=cvvn)+
  stat_ecdf()+
  #facet_wrap(~variable)+
  #geom_vline(xintercept = c(0.03,0.01),color="grey40",linetype=3)+
  theme_classic()+
  scale_y_continuous(name="ECDF")+
  scale_x_continuous(name="vote",trans="log10")+
  scale_linetype_manual(name="",values=c(1,3))+
  annotate("rect", xmin = 0.003, xmax = 0.01, ymin = 0, ymax = 1, alpha = .2)+
  scale_color_brewer(palette = "Dark2",name="")
ggsave("votes.pdf",width=7,height=5.5)


ggplot(aes(x=vote,color=variable,linetype=value),data=cvvp)+
  stat_bin(data=subset(cvvp,variable=="TP"),aes(y=cumsum(after_stat(count))),geom="step")+
  stat_bin(data=subset(cvvp,variable=="FP"),aes(y=cumsum(after_stat(count))),geom="step")+
  theme_classic()+
  scale_y_continuous(name="ECDF")+
  scale_x_continuous(name="vote",trans="log10")+
  scale_linetype_manual(name="",values=c(1,3))+
  annotate("rect", xmin = 0.003, xmax = 0.01, ymin = 0, ymax = 1, alpha = .1)+
  annotate("rect", xmin = 0.003, xmax = 0.03, ymin = 0, ymax = 1, alpha = .1)+
  scale_color_brewer(palette = "Dark2",name="")
ggsave("votes.pdf",width=7,height=5.5)


ggplot(aes(x=log10(vote),y=vote_n),data=cvvp)+
  facet_grid(value~variable)+
  geom_bin_2d()+
  theme_classic()+
  scale_y_continuous(name="normalized vote")+
  coord_cartesian(ylim= c(0.4995,1.001))+
  #geom_vline(xintercept = 0.03,color="grey40")+
  scale_x_continuous(name="vote")+
  scale_fill_continuous(trans="log10")
ggsave("votes-norm.pdf",width=7,height=7)


ggplot(aes(x=vote_n,color=variable,linetype=value),data=cvvp)+
  stat_ecdf()+facet_wrap(~value)+
  theme_classic()+
  scale_y_continuous(name="ECDF")+
  scale_x_continuous(name="vote (normalized)")+
  scale_linetype_manual(name="",values=c(1,3))+
  scale_color_brewer(palette = "Dark2",name="")+
  coord_cartesian(xlim= c(0.4995,1.001))+
  theme(legend.position = "none")
ggsave("votesnorm.pdf",width=8.8,height=5.5)

  e = 0.3 * diff(range(cvv$vote))
  
  dens = density(cvv$vote, adjust=0.2, from=min(cvv$vote)-e, to=cvv$vote +e)
  dens = data.frame(x=cvv$vote, y=cvv$vote)
  
  # Plot kernel density (blue), ecdf (red) and smoothed ecdf (black)
  ggplot(aes(x=vote, y=cumsum(y)/sum(y))) + 
    geom_line( size=0.7,) +
    stat_ecdf(colour="red", size=0.6, alpha=0.6) +
    theme_classic() +
    labs(title=paste0("adj=",adj))



