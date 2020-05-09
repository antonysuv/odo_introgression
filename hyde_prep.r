library('ape')
library('phytools')
library('phangorn')
library('Rtsne')
library('MASS')
library('gridExtra')
library('ggplot2')
library('pals')
library('reshape2')



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
#Aeshnoidea=tt[50:63]
Cordulegastroidea=tt[64:69]
Libelluloidea=tt[70:83]
Platystictidae=tt[43]
Aeshnidae=tt[50:56]
Gomphidae_Petaluridae=tt[57:63]
Libellulidae=tt[75:83]


total=read.table("hyde_all_tri-out_noephemera.txt",header=T)
total=total[complete.cases(total),]
total$Pvalue=p.adjust(total$Pvalue,method="bonferroni")
total$D=(total$ABBA-total$ABAB)/(total$ABBA+total$ABAB)
total=total[total$Gamma<=1 & total$Gamma>=0, ]

#Chisq test for D significance
chiPd=unlist(lapply(apply(total[,c("ABBA","ABAB")],1,chisq.test),"[[","p.value"))
chiPd=p.adjust(chiPd,method="bonferroni")
total$PvalueD=chiPd



total$P1subo=ifelse(total$P1 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P1 %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$Hybridsubo=ifelse(total$Hybrid %in% Zygoptera,"Zygoptera",
                    ifelse(total$Hybrid %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$P2subo=ifelse(total$P2 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P2 %in% Anisoptera,"Anisoptera","Anisozygoptera"))
    

total$P1sub=ifelse(total$P1 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P1 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P1 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P1 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P1 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P1 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P1 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P1 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P1 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))

total$Hybridsub=ifelse(total$Hybrid %in% Lestoidea,"Lestoidea",
                   ifelse(total$Hybrid %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$Hybrid %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$Hybrid %in% Platystictidae,"Platystictidae",
                   ifelse(total$Hybrid %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$Hybrid %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$Hybrid %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$Hybrid %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$Hybrid %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))


total$P2sub=ifelse(total$P2 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P2 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P2 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P2 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P2 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P2 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P2 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P2 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P2 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))



nodups=!apply(t(apply(total[,c("P1sub","Hybridsub","P2sub")],1,duplicated)),1,any)


#total$Order=ifelse((total$P2subo=="Zygoptera" & total$Hybridsubo=="Zygoptera" &  (total$P1subo=="Anisoptera" | total$P1subo=="Anisozygoptera" )) | ( (total$P2subo=="Anisoptera" | total$P2subo=="Anisozygoptera") & (total$Hybridsubo=="Anisoptera" total$P1subo=="Anisozygoptera") &  total$P1subo=="Zygoptera"),"Odonata","RANDOM")
#total$Total="All"

total$Order="Odonata"


total$Suborder=ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Zygoptera),1,all) & nodups,"Zygoptera",
                      ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Anisoptera),1,all) & nodups,"Anisoptera",
                     ifelse(total$Hybrid==Anisozygoptera & nodups & total$P1subo!=total$P2subo,"Anisozygoptera","RANDOM")))

total$Superfamily=ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM")))))       
                              


total$focalclade=ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Aeshnidae),1,all),"Aeshnidae",
                      ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Gomphidae_Petaluridae),1,all),"Gomphidae+Petaluridae",
                      ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Libellulidae),1,all),"Libellulidae","RANDOM")))

total=total[total$D >= 0,]
total$introgressionid=ifelse(total$PvalueD<0.05 & total$Pvalue<10^-6,"Introgression","None") 
sunset=c("#364B9A" ,"#4A7BB7", "#6EA6CD", "#98CAE1" ,"#C2E4EF" ,"#EAECCC", "#FEDA8B" ,"#FDB366", "#F67E4B", "#DD3D2D", "#A50026")
#D and Gamma distributions Violin plots

total_ord=melt(total[,c("D","Gamma","Pvalue","Order","Suborder","Superfamily","focalclade","PvalueD")],value.name="taxon",id=c("D","Gamma","Pvalue","PvalueD"))
total_ord=total_ord[total_ord$taxon!="RANDOM",]
total_ord$variable=replace(as.character(total_ord$variable),as.character(total_ord$variable)=="focalclade","Focal clade")
total_ord$variable_f=factor(total_ord$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))

a1=ggplot(total_ord[total_ord$PvalueD<0.05 & total_ord$D > 0,], aes(x=taxon, y=D))+geom_violin(fill='olivedrab3')+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1),plot.margin=unit(c(0,0.1,0,0.1), "cm"))+ylab("D")+xlab("")+ggtitle("A")

a2=ggplot(total_ord[total_ord$Pvalue<10^-6,], aes(x=taxon, y=Gamma))+geom_violin(fill='cornflowerblue')+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1),plot.margin=unit(c(0,0.1,0,0.1), "cm"))+ylab(expression(gamma))+xlab("")+ggtitle("B")

dat=total_ord[total_ord$PvalueD<0.05 & total_ord$D > 0 & total_ord$Pvalue<10^-6,]
dat$density = get_density(dat$Gamma, dat$D, n = 100)

a4=ggplot(dat,aes(Gamma, D, color = density)) + geom_point(size=0.4,stroke=0)+geom_smooth(method = "auto", size = 0.5,color="black")+scale_color_gradientn(colors = sunset[1:6],guide=guide_colorbar(barwidth = 6,barheight =0.5,label.theme = element_text(size=8)),name="")+labs(x = expression(gamma),y="D")+theme(legend.position=c(0.5,1.05) ,legend.direction="horizontal",legend.background = element_blank())+ggtitle("D")


total_p=melt(total[,c("introgressionid","Order","Suborder","Superfamily","focalclade")],id="introgressionid",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="focalclade","Focal clade")
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))

a3=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=introgressionid))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5))+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1),legend.position=c(0.5,1.4) ,legend.direction="horizontal",legend.background = element_blank())+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("wheat","grey50"),name="")+ggtitle("C")




#tSNE
totalsig=total[total$introgressionid=="Introgression",c("Gamma","D")]
total_dr=total[total$PvalueD<0.05 & total$D > 0 & total$Pvalue<10^-6,c(7:21)]/apply(total[total$PvalueD<0.05 & total$D > 0 & total$Pvalue<10^-6,c(7:21)],1,sum)
tsne_out=Rtsne(total_dr)
dat1=data.frame(tsne_out$Y)
names(dat1)=c("tSNE1","tSNE2")
#dat1=cbind(dat1,introgressionid=ifelse(total$Pvalue<10^-6,"Introgression","None"))
#dat1=cbind(dat1,introgressionidD=ifelse(total$PvalueD<0.05,"Introgression","None"))
dat1=cbind(dat1,totalsig)


a5=ggplot(dat1) + geom_point(aes(tSNE1, tSNE2, color = D),size=0.4,stroke=0)+scale_color_gradientn(colors = sunset,name="D",guide=guide_colorbar(barwidth = 6,barheight =0.5,label.theme = element_text(size=8)))+theme_classic()+theme(legend.position=c(0.5,1.05) ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank())+ggtitle("E")

a6=ggplot(dat1) + geom_point(aes(tSNE1, tSNE2, color = Gamma),size=0.3,stroke=0)+scale_color_gradientn(colors = sunset,name=expression(gamma),guide=guide_colorbar(barwidth = 6,barheight =0.5,label.theme = element_text(size=8)))+theme_classic()+theme(legend.position=c(0.5,1.05) ,legend.direction="horizontal",legend.background = element_blank())+ggtitle("F")




lay = rbind(c(1,1,1),
             c(2,2,2),
           c(3,3,3),
            c(4,5,6))

quartz(width=7,height=9.5) 
grid.arrange(a1,a2,a3,a4,a5,a6,layout_matrix = lay)

quartz.save("hyde_distribution.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("hyde_distribution.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)





#Comparisons
#D means
total_Dsig=total[total$PvalueD<0.05,]
total_Gsig=total[total$Pvalue<10^-6,]

get_wilx=function(m,stats)
{
    m_out=c()
    for(i in c("Order","Suborder","Superfamily","focalclade"))
    {
       
        for (cl in unique(m[,i]))
        {
            
            my_mean=mean(m[m[,i]==cl,stats])
            my_pg=wilcox.test(m[m[,i]==cl,stats],m[,stats],alternative="greater")$p.value
            my_pt=wilcox.test(m[m[,i]==cl,stats],m[,stats],alternative="two.sided")$p.value
            my_pl=wilcox.test(m[m[,i]==cl,stats],m[,stats],alternative="less")$p.value
            my_t=t.test(m[m[,i]==cl,stats],mu=0.5)$p.value
            m_out=rbind(m_out,c(i,cl,my_mean,my_pg,my_pt,my_pl,my_t))
        }
            
    }
    m_out=as.data.frame(m_out)
    names(m_out)=c("rank","taxon","average","p_greater","p_two","p_less","p_ttest")
    
    return(m_out)
        
}

#Older quartets have smaller D

young_t=total_Dsig[total_Dsig$Superfamily=="Coenagrionoidea"| total_Dsig$focalclade=="Libellulidae" | total_Dsig$focalclade=="Aeshnidae" ,"D"]
old_t=total_Dsig[total_Dsig$P1subo=="Zygoptera" & total_Dsig$P2subo=="Anisoptera" | (total_Dsig$P1subo=="Anisoptera" & total_Dsig$P2subo=="Zygoptera"),"D"]
wilcox.test(young_t,old_t,alternative="greater")


















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
