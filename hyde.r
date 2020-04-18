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


################################################################ RUN HyDe ######################################################################## 



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
total=total[complete.cases(total),]
total$Pvalue=p.adjust(total$Pvalue,method="bonferroni")
total$D=(total$ABBA-total$ABAB)/(total$ABBA+total$ABAB)
total=total[total$Gamma<=1 & total$Gamma>=0, ]

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
                              

#Dimentionality reduction
total_dr=total[,c(7:21)]/apply(total[,c(7:21)],1,sum)
phate_out=phate(total_dr)
umap_out=umap(total_dr)
tsne_out=Rtsne(total_dr)

#Main dr figure
dat1=data.frame(tsne_out$Y)
names(dat1)=c("tSNE1","tSNE2")
introg_d=cbind(dat1,introgressionid=ifelse(total$Pvalue<10^-6,"Introgression","None"))
dat1=cbind(dat1,total[,c("Gamma","D")])
t1=ggplot(introg_d) + geom_point(aes(tSNE1, tSNE2, color = introgressionid),size=0.001)+ scale_color_manual(values=c("#5cb85c","#337ab7"),name="")+theme_classic()+ guides(colour = guide_legend(override.aes = list(size=3)))+ theme(legend.position=c(0.15, 0.1),legend.background = element_blank())
t2=ggplot(dat1) + geom_point(aes(tSNE1, tSNE2, color = D),size=0.01)+scale_color_gradientn(colors = tol(12),name="")+theme_classic()+theme(legend.position=c("bottom") ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+ ggtitle("D")
t3=ggplot(dat1[total$Pvalue<10^-6,]) + geom_point(aes(tSNE1, tSNE2, color = Gamma),size=0.01)+scale_color_gradientn(colors = tol(12),name="")+theme_classic()+theme(legend.position=c("bottom") ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle(expression(gamma))
lay = rbind(c(1,1,1,NA),
             c(2,2,3,3))
             
quartz(width=5.9,height=8.4) 
grid.arrange(t1,t2,t3,layout_matrix = lay)
quartz.save("hyde_tsne.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("hyde_tsne.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)

#Suppl dr figure 
#PHATE
dat1=data.frame(phate_out$embedding)
names(dat1)=c("PHATE1","PHATE2")
introg_d=cbind(dat1,introgressionid=ifelse(total$Pvalue<10^-6,"Introgression","None"))
dat1=cbind(dat1,total[,c("Gamma","D")])
t4=ggplot(introg_d) + geom_point(aes(PHATE1,PHATE2, color = introgressionid),size=0.001)+ scale_color_manual(values=c("#5cb85c","#337ab7"),name="")+theme_classic()+ guides(colour = guide_legend(override.aes = list(size=3)))+ theme(legend.position=c(0.15, 0.1),legend.background = element_blank())+xlim(-0.01,0.1)+xlim(-0.03,0.04)
t5=ggplot(dat1) + geom_point(aes(PHATE1,PHATE2, color = D),size=0.01)+scale_color_gradientn(colors = tol(12),name="")+theme_classic()+theme(legend.position=c("bottom") ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+ ggtitle("D")+xlim(-0.03,0.04)
t6=ggplot(dat1[total$Pvalue<10^-6,]) + geom_point(aes(PHATE1,PHATE2, color = Gamma),size=0.01)+scale_color_gradientn(colors = tol(12),name="")+theme_classic()+theme(legend.position=c("bottom") ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle(expression(gamma))+xlim(-0.03,0.04)

#UMAP
dat1=data.frame(umap_out$layout)
names(dat1)=c("UMAP1","UMAP2")
introg_d=cbind(dat1,introgressionid=ifelse(total$Pvalue<10^-6,"Introgression","None"))
dat1=cbind(dat1,total[,c("Gamma","D")])
t7=ggplot(introg_d) + geom_point(aes(UMAP1,UMAP2, color = introgressionid),size=0.001)+ scale_color_manual(values=c("#5cb85c","#337ab7"),name="")+theme_classic()+ guides(colour = guide_legend(override.aes = list(size=3)))+ theme(legend.position=c(0.15, 0.1),legend.background = element_blank())
t8=ggplot(dat1) + geom_point(aes(UMAP1,UMAP2, color = D),size=0.01)+scale_color_gradientn(colors = tol(12),name="")+theme_classic()+theme(legend.position=c("bottom") ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+ ggtitle("D")
t9=ggplot(dat1[total$Pvalue<10^-6,]) + geom_point(aes(UMAP1,UMAP2, color = Gamma),size=0.01)+scale_color_gradientn(colors = tol(12),name="")+theme_classic()+theme(legend.position=c("bottom") ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle(expression(gamma))

lay = rbind(c(1,1,1,NA,4,4,4,NA),
             c(2,2,3,3,5,5,6,6))

quartz(width=10.7,height=9.5) 
grid.arrange(t4,t5,t6,t7,t8,t9,layout_matrix = lay)
quartz.save("hyde_phate_umap_suppl.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("hyde_phate_umap_suppl.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)


#STATS for Focal clades 

#Epiophlebia 
total=total[total$Pvalue<10^-6,]
total[total$Hybrid=="P49_Epiophlebia_superstes" & total$P1 %in% Lestoidea  &  total$P2 %in% Aeshnoidea[1:7],]$Gamma


#D and Gamma distributions Violin plots

total_ord=melt(total[,c("D","Gamma","Order")])
a1=ggplot(total_ord[total_ord$Order!="RANDOM",], aes(x=Order, y=value,fill=variable))+geom_violin()+labs(x="Order", y = "Value")+scale_fill_manual(values=c("Grey", "goldenrod2"),name="",labels=c("D",expression(gamma)))+facet_wrap(~variable)+ theme(axis.text.x = element_text(size = 8))+ geom_boxplot(width=0.01,outlier.size=-1)+ggtitle("A")+stat_summary(fun.y=median, geom="point", size=2, color="black")

total_suborder=melt(total[,c("D","Gamma","Suborder")])
a2=ggplot(total_suborder[total_suborder$Suborder!="RANDOM",], aes(x=Suborder, y=value,fill=variable))+geom_violin()+labs(x="Suborder", y = "Value")+scale_fill_manual(values=c("Grey", "goldenrod2"),name="",labels=c("D",expression(gamma)))+facet_wrap(~variable)+ theme(axis.text.x = element_text(size = 8))+ geom_boxplot(width=0.01,outlier.size=-1)+ggtitle("B")+stat_summary(fun.y=median, geom="point", size=2, color="black")

total_superfam=melt(total[,c("D","Gamma","Superfamily")])
a3=ggplot(total_superfam[total_superfam$Superfamily!="RANDOM",], aes(x=Superfamily, y=value,fill=variable))+geom_violin()+labs(x="Superfamily", y = "Value")+scale_fill_manual(values=c("Grey", "goldenrod2"),name="",labels=c("D",expression(gamma)))+facet_wrap(~variable)+ theme(axis.text.x = element_text(size = 8,angle=10))+ geom_boxplot(width=0.01,outlier.size=-1)+ggtitle("C")+stat_summary(fun.y=median, geom="point", size=2, color="black")

grid.arrange(a1,a2,a3,ncol=1,nrow=3)






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

