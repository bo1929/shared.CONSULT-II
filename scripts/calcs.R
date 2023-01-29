require(ggplot2)
require(scales)
SL = 32
p = 3

f = function(x,k,l,SL) { 1 - (1 - (1 - x/SL)^k)^l }

f(3,15,2,32)

K = function(l,p,alpha,L) {round(log(1-(1-alpha)^(1/l))/log(1-p/L));}

K(4,2,0.95,10)

alpha=0.95
ggplot(data.frame(x = c(1, 3*p)), aes(x)) + theme_classic()+
  mapply(function( L) {
    stat_function(fun = f, args = list(k = K(l=L,p=p,alpha=alpha,L=SL), l=L, SL=SL), aes(color=paste(L,K(l=L,p=p,alpha=alpha,L=SL),sep="/")))
  },L=c(7:14))+
  geom_vline(xintercept = p,color="black",linetype=2)+
  geom_hline(yintercept = alpha,color="black",linetype=2)+
  scale_color_brewer(palette = "Spectral",name="l/h")+ylab("probability of a match")+xlab("Hamming distance")+scale_linetype(name="l")
ggsave("clus-95.pdf",width=5,height = 4)

alpha=0.48
ggplot(data.frame(x = c(1, 3*p)), aes(x)) + theme_classic()+
  mapply(function( L) {
    stat_function(fun = f, args = list(k = K(l=L,p=p,alpha=alpha,L=SL), l=L, SL=SL), aes(color=paste(L,K(l=L,p=p,alpha=alpha,L=SL),sep="/")))
  },L=c(1:8))+
  geom_vline(xintercept = p,color="black",linetype=2)+
  geom_hline(yintercept = alpha,color="black",linetype=2)+
  scale_color_brewer(palette = "Spectral",name="l/h")+ylab("probability of a match")+xlab("Hamming distance")+scale_linetype(name="l")
ggsave("clus-10.pdf",width=5,height = 4)

alpha=0.1
ggplot(data.frame(x = c(1, 3*p)), aes(x)) + theme_classic()+
  mapply(function( L) {
    stat_function(fun = function(x,k,l,SL) (150-SL+1)*f(x,k,l,SL), args = list(k = K(l=L,p=p,alpha=alpha,L=SL), l=L, SL=SL), aes(color=paste(L,K(l=L,p=p,alpha=alpha,L=SL),sep="/")))
  },L=c(2:8))+
  geom_vline(xintercept = 4,color="black",linetype=2)+
  geom_hline(yintercept = 3,color="black",linetype=2)+
  scale_color_brewer(palette = "Dark2",name="l/k")+ylab("probability of a match")+xlab("distance")+scale_linetype(name="l")



alpha=0.1
ggplot(data.frame(x = c(1, 16)), aes(x)) + theme_classic()+
  mapply(FUN=function( L, K) {
    stat_function(fun = function(x,k,l,SL) (1+150-SL)*f(x,k,l,SL), args = list(k = K, l=L, SL=SL), aes(color=as.factor(K),linetype=as.factor(L)))
  },L=rep(c(2,4,1),each=5), K=c(6:10)*2-1)+
  #geom_vline(xintercept = 3,color="black",linetype=2)+
  #geom_vline(xintercept = SL*0.4,color="black",linetype=2)+
  #geom_hline(yintercept = 90,color="black",linetype=2)+
  scale_y_log10(lim=c(0.1,100))+
  scale_linetype_manual(values=c(2,1,3),name=expression(l))+
  scale_x_continuous(labels=function(x) paste(x,percent(x/SL),sep="\n"))+
  scale_color_brewer(palette = "Dark2",name="h")+ylab("Expected number of matches")+xlab("Hamming distance")+
  stat_function(fun = function(x) (1+150-SL)*f(x,35-7,1,35)/2, color="black",aes(linetype=as.factor(1)))
ggsave("exp.pdf",width=5,height = 4)

alpha=0.1
ggplot(data.frame(x = c(1, 16)), aes(x)) + theme_classic()+
  mapply(FUN=function( L, K) {
    stat_function(fun = function(x,k,l,SL) (1+150-SL)*f(x,k,l,SL), 
                  args = list(k = K, l=L, SL=SL), 
                  aes(color=as.factor(K),linetype=as.factor(L)))
  },L=rep(c(2),each=5), K=c(6,10,14))+
  #geom_vline(xintercept = 3,color="black",linetype=2)+
  #geom_vline(xintercept = SL*0.4,color="black",linetype=2)+
  #geom_hline(yintercept = 90,color="black",linetype=2)+
  scale_y_continuous(lim=c(0.1,100))+
  scale_linetype_manual(values=c(1,2,3),name=expression(h))+
  scale_x_continuous(labels=function(x) paste(x,percent(x/SL),sep="\n"))+
  scale_color_manual(name="l",values = c("gray50","black","gray80"))+
  ylab("Expected number of kmer matches per 150bp read")+
  xlab("Hamming distance")+
  #stat_function(fun = function(x) (1+150-SL)*f(x,35-7,1,35)/2, color="blue",linetype=1)+
  geom_vline(xintercept = 3,color="red",linetype=3)+
  geom_vline(xintercept = 10,color="red",linetype=3)+
  theme(legend.position = "none")
  ggsave("exp-h.pdf",width=5,height = 4)


hueh = function(Mem,N) log2(Mem/(8*N))-1

hueh = function(N,b=7) ceiling(1/2 * log2(7/8*N/b))
frl = function(x,k,SL,p) { round((log2(1 - p)/log2(1 - (1 - x/SL)^k) ))}
hueh(2^33)
frl(3,15,32,0.45)




ggplot(data.frame(x = c(1, 3*p)), aes(x)) + theme_classic()+
  mapply(function(N) {
    stat_function(fun = f, args = list(k = hueh(N), l=frl(3,hueh(N),32,0.5), SL=SL), 
                  aes(color=paste(frl(3,hueh(N),32,0.5),hueh(N),N,sep="/")))
  },N=c(2^22,2^28,2^30,2^33))+
  geom_vline(xintercept = 3,color="black",linetype=2)+
  geom_hline(yintercept = 0.5,color="black",linetype=2)+
  scale_color_brewer(palette = "Spectral",name="l/h")+
  ylab("probability of a match")+xlab("Hamming distance")+scale_linetype(name="l")


ggplot(data.frame(x = c(2^20, 2^33.1)), aes(x)) + theme_classic()+
  mapply(function(b) {
    stat_function(fun = hueh, args = list(b=b), 
          aes(color=b,linetype="h"))
  },b=c(4:20))+
  mapply(function(b) {
    stat_function(fun = function(x,b) frl(3,hueh(x,b=b),32,0.45), args = list(b=b), 
                  aes(color=b,linetype="l"))
  },b=c(4:20))+
  #scale_color_brewer(palette = "Set2",name="b")+
  geom_vline(xintercept = 2^33,color="grey40",linetype=2)+
  scale_y_continuous("setting",breaks=c(2,5,8,11,14))+
  scale_x_continuous(name="k-mer count (G)",trans="log10")+
  scale_color_viridis_c(guide = guide_legend(breaks=2),breaks=(2:10)*2)+
  theme(legend.position = c(0.385,.92),legend.direction = "horizontal",
        legend.margin = margin(0),legend.box.margin = margin(0),
        legend.spacing = unit(0,"pt")
  )+
  scale_linetype(name="")
ggsave("selected-h-l-perb.pdf",width=4.2,height = 3.7)

gl = function (N,b) 4 * b* 2^(2 * hueh(N,b=b)) * frl(3,hueh(N,b=b),32,0.45) + 8 * N
ggplot(data.frame(x = c(2^20, 2^33.2)), aes(x)) + theme_classic()+
  mapply(function(b) { stat_function(fun = function(x) gl(x,b=b),aes(color=b) )},
         b=c(4:20))+
  #scale_color_brewer(palette = "Spectral",name="b")+
  geom_vline(xintercept = 2^33,color="grey40",linetype=2)+
  geom_hline(yintercept = 2^37,color="grey40",linetype=2)+
  scale_y_continuous("Memory (GB)",trans="log10",breaks=2^(c(11:19)*2+1),labels = number_bytes)+
  scale_x_continuous(name="k-mer count (G)",trans="log10")+
  scale_color_viridis_c(guide = guide_legend(breaks=2),breaks=(2:10)*2)+
  theme(legend.position = c(0.44,.95),legend.direction = "horizontal",
        legend.margin = margin(0),legend.box.margin = margin(0),
        legend.spacing = unit(0,"pt")
        )
ggsave("mem-perb.pdf",width=4.2,height = 3.7)
