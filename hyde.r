library('ape')
library('phytools')
library('phangorn')
library('Rtsne')
library('MASS')
library('gridExtra')
library('ggplot2')
library('pals')
library('reshape2')
library('umap')
library("phateR")


get_density = function(x, y, ...) 
{
  dens = MASS::kde2d(x, y, ...)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}


list_triples=function(taxa)
{
    res=rbind()
    ks = taxa
    for (i in 1:(length(ks)-2))
    {
        for (j in (i+1):(length(ks)-1))
        {
           for (k in (j+1):length(ks))
           {
               res=rbind(res,c(ks[i],ks[j], ks[k]))
               res=rbind(res,c(ks[i],ks[k], ks[j]))
               res=rbind(res,c(ks[j],ks[i], ks[k]))
               
            }
        }    
    }
    write.table(res,"tri.txt",quote=F,row.names=F,col.names=F)
    
}        
    
taxa_map=function(tree)
{
    phy=read.tree(tree)
    d=data.frame(phy$tip.label,paste("P",1:length(phy$tip.label),"_",phy$tip.label,sep=""))
    write.table(d,"taxa_map.txt",quote=F,row.names=F,col.names=F)
}
    
tripls_from_taxa=function(taxa,taxa_map)
{
    d=read.table(taxa_map)
    taxa=as.vector(d[d$V1 %in% taxa,"V2"])
    list_triples(taxa)
}


taxa_map("BUSCO50_dna_pasta_nopart_iqtree_root.tre")
tripls_from_taxa(gg,"taxa_map.txt") 



tt=read.table("tri.txt")





#Violin plots
tt=as.character(read.table("taxa_map.txt")$V2)
Zygoptera=tt[1:48]
Anisoptera=tt[50:83]
Anisozygoptera=tt[49]
Lestoidea=tt[44:48]
Calopterygoidea=tt[1:19]
Coenagrionoidea=tt[20:42]
Aeshnoidea=tt[50:63]
Cordulegastroidea=tt[64:69]
Libelluloidea=tt[70:83]

total=read.table("hyde_all_tri-out_noephemera.txt",header=T)
total$Pvalue=p.adjust(total$Pvalue,method="bonferroni")
total=total[total$Pvalue<10^-6,]
total$D=(total$ABBA-total$ABAB)/(total$ABBA+total$ABAB)
total$Order="Odonata"

total$Suborder=ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Zygoptera),1,all),"Zygoptera",
                      ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Anisoptera),1,all),"Anisoptera",
                      ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Anisozygoptera),1,any),"Anisozygoptera","RANDOM")))

total$Superfamily=ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Aeshnoidea),1,all),"Aeshnoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM"))))))       
                              




total_ord=melt(total[,c("D","Gamma","Order")])
a1=ggplot(total_ord[total_ord$Order!="RANDOM",], aes(x=Order, y=value,fill=variable))+geom_violin()+labs(x="Order", y = "Value")+scale_fill_manual(values=c("Grey", "goldenrod2"),name="",labels=c("D",expression(gamma)))+facet_wrap(~variable)+ theme(axis.text.x = element_text(size = 8))+ geom_boxplot(width=0.01,outlier.size=-1)+ggtitle("A")+stat_summary(fun.y=median, geom="point", size=2, color="black")

total_suborder=melt(total[,c("D","Gamma","Suborder")])
a2=ggplot(total_suborder[total_suborder$Suborder!="RANDOM",], aes(x=Suborder, y=value,fill=variable))+geom_violin()+labs(x="Suborder", y = "Value")+scale_fill_manual(values=c("Grey", "goldenrod2"),name="",labels=c("D",expression(gamma)))+facet_wrap(~variable)+ theme(axis.text.x = element_text(size = 8))+ geom_boxplot(width=0.01,outlier.size=-1)+ggtitle("B")+stat_summary(fun.y=median, geom="point", size=2, color="black")

total_superfam=melt(total[,c("D","Gamma","Superfamily")])
a3=ggplot(total_superfam[total_superfam$Superfamily!="RANDOM",], aes(x=Superfamily, y=value,fill=variable))+geom_violin()+labs(x="Superfamily", y = "Value")+scale_fill_manual(values=c("Grey", "goldenrod2"),name="",labels=c("D",expression(gamma)))+facet_wrap(~variable)+ theme(axis.text.x = element_text(size = 8,angle=10))+ geom_boxplot(width=0.01,outlier.size=-1)+ggtitle("C")+stat_summary(fun.y=median, geom="point", size=2, color="black")

grid.arrange(a1,a2,a3,ncol=1,nrow=3)


#Dimentionality reduction
total_dr=total[,c(7:21)]/apply(total[,c(7:21)],1,sum)
phate_out=phate(total_dr)
umap_out=umap(total_dr)
tsne_out=Rtsne(total_dr)

dat1=data.frame(phate_out$embedding)
dat1$Gamma=total$Gamma
dat1$D=total$D

dat2=data.frame(umap_out$layout)
names(dat2)=c("UMAP1","UMAP2")
dat2$Gamma=total$Gamma
dat2$D=total$D

dat3=data.frame(tsne_out$Y)
names(dat3)=c("tSNE1","tSNE2")
dat3$Gamma=total$Gamma
dat3$D=total$D

g1=ggplot(dat1) + geom_point(aes(PHATE1, PHATE2, color = Gamma),size=0.01)+scale_color_gradientn(colors = viridis(30),name=expression(gamma))+theme_classic()+ggtitle("A")
g2=ggplot(dat1) + geom_point(aes(PHATE1, PHATE2, color = D),size=0.01)+scale_color_gradientn(colors = viridis(30))+theme_classic()+ggtitle("")
g3=ggplot(dat2) + geom_point(aes(UMAP1, UMAP2, color = Gamma),size=0.01)+scale_color_gradientn(colors = viridis(30),name=expression(gamma))+theme_classic()+ggtitle("B")
g4=ggplot(dat2) + geom_point(aes(UMAP1, UMAP2, color = D),size=0.01)+scale_color_gradientn(colors = viridis(30))+theme_classic()+ggtitle("")
g5=ggplot(dat3) + geom_point(aes(tSNE1, tSNE2, color = Gamma),size=0.01)+scale_color_gradientn(colors = viridis(30),name=expression(gamma))+theme_classic()+ggtitle("C")
g6=ggplot(dat3) + geom_point(aes(tSNE1, tSNE2, color = D),size=0.01)+scale_color_gradientn(colors = viridis(30))+theme_classic()+ggtitle("")

grid.arrange(g1,g2,g3,g4,g5,g6,ncol=2,nrow=3)
quartz.save("hyde_phate_total_sitespatterns", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)

pca_out=prcomp(total_dr)
dat4=data.frame(pca_out$x[,1:2])
dat4$Gamma=total$Gamma
dat4$D=total$D
ggplot(dat4) + geom_point(aes(PC1,PC2, color = Gamma),size=0.01)+scale_color_gradientn(colors = viridis(30),name=expression(gamma))+theme_classic()
ggplot(dat4) + geom_point(aes(PC1,PC2, color = D),size=0.01)+scale_color_gradientn(colors = viridis(30))+theme_classic()




#Total tree test
total=read.table("hyde_all_tri-out_noephemera.txt",header=T)
total=total[complete.cases(total), ]
total$D=(total$ABBA-total$ABAB)/(total$ABBA+total$ABAB)
total=total[total$Gamma > 0 & total$Gamma < 1, ]
total$Pvalue=p.adjust(total$Pvalue,method="bonferroni")

dat=total[,c("Gamma","D","Pvalue")]
g1=ggplot(dat)+geom_histogram(aes(Gamma),bins=500,col="black")+labs(x = expression(gamma),y="N quartets")+ggtitle("A")
g2=ggplot(dat)+geom_histogram(aes(D),bins=500,col="black")+labs(x = "D",y="N quartets")+ggtitle("")
dat$density = get_density(dat$Gamma, dat$D, n = 100)
g3=ggplot(dat) + geom_point(aes(Gamma, D, color = density),size=0.1)+scale_color_gradientn(colors = viridis(30))+labs(x = expression(gamma),y="D")+ggtitle("")
g4=ggplot(dat) + geom_point(aes(Gamma, D, color = Pvalue),size=0.1)+scale_color_gradientn(colors = viridis(30))+labs(x = expression(gamma),y="D")+ggtitle("")

dat=dat[dat$Pvalue<10^-6,c("Gamma","D","Pvalue")]
g5=ggplot(dat)+geom_histogram(aes(Gamma),bins=500,col="black")+labs(x = expression(gamma),y="N quartets")+ggtitle("B")
g6=ggplot(dat)+geom_histogram(aes(D),bins=500,col="black")+labs(x = "D",y="N quartets")+ggtitle("")
dat$density = get_density(dat$Gamma, dat$D, n = 100)
g7=ggplot(dat) + geom_point(aes(Gamma, D, color = density),size=0.1)+scale_color_gradientn(colors = viridis(30))+labs(x = expression(gamma),y="D")+ggtitle("")
g8=ggplot(dat) + geom_point(aes(Gamma, D, color = Pvalue),size=0.1)+scale_color_gradientn(colors = viridis(30))+labs(x = expression(gamma),y="D")+ggtitle("")

grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=4,nrow=2)
quartz.save("hyde_total", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
