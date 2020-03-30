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






#Epiophlebia test 
Aeshnidae=c("P50_Aeshna_palmata","P51_Anax_junius","P52_Anax_walsinghami","P53_Anax_parthenope","P54_Gynacantha_tibiata","P55_Austroaeschna_subapicalis","P56_Telephlebia_godeffroy")
Lestoidea=c("P43_Protosticta_beaumonti","P44_Archilestes_grandis","P45_Indolestes_peregrinus","P46_Episynlestes_cristatus","P47_Synlestes_weyersii","P48_Perissolestes_remotus")
Outgroup=c("P84_Ephemera_danica","P85_Isonychia_kiangsinensis")
write.table(tt[tt$V2=="P49_Epiophlebia_superstes" & tt$V1 %in% Lestoidea & tt$V3 %in% Aeshnidae,],"epio_tri.txt",col.names=F,row.names=F,quote=F)
write.table(tt[tt$V2=="P49_Epiophlebia_superstes" & !tt$V1 %in% Lestoidea & !tt$V3 %in% Aeshnidae & !tt$V1 %in% Outgroup & !tt$V3 %in% Outgroup,],"epio_tri_random.txt",col.names=F,row.names=F,quote=F)



epio_hyde=read.table("hyde_epio_tri-out.txt",header=T)
epio_hyde$Pvalue=p.adjust(epio_hyde$Pvalue,method="bonferroni")
epio_hyde=epio_hyde[epio_hyde$Pvalue<0.05,]
epio_hyde$Gamma[which(epio_hyde$Gamma > 1)] = 1
epio_hyde$Gamma[which(epio_hyde$Gamma < 0)] = 0
epio_hyde$D=(epio_hyde$ABBA-epio_hyde$ABAB)/(epio_hyde$ABBA+epio_hyde$ABAB)


epio_hyde_random=read.table("hyde_epio_tri_random-out.txt",header=T)
epio_hyde_random$Pvalue=p.adjust(epio_hyde_random$Pvalue,method="bonferroni")
epio_hyde_random=epio_hyde_random[epio_hyde_random$Gamma > 0 & epio_hyde_random$Gamma < 1 & epio_hyde_random$Pvalue < 0.05, ]
epio_hyde_random$D=(epio_hyde_random$ABBA-epio_hyde_random$ABAB)/(epio_hyde_random$ABBA+epio_hyde_random$ABAB)

total=read.table("hyde_all_tri-out_noephemera.txt",header=T)
total$Pvalue=p.adjust(total$Pvalue,method="bonferroni")
total=total[total$Gamma > 0 & total$Gamma < 1 & total$Pvalue<0.05, ]
total$D=(total$ABBA-total$ABAB)/(total$ABBA+total$ABAB)
#Violin plot
vp_hydeG=melt(list(Epiophlebia=epio_hyde$Gamma,Random=epio_hyde_random$Gamma,Total=total$Gamma))
vp_hydeG$stat="Gamma"
vp_hydeD=melt(list(Epiophlebia=epio_hyde$D,Random=epio_hyde_random$D,Total=total$D))
vp_hydeD$stat="D"
vp_hyde1=rbind(vp_hydeG,vp_hydeD)
ggplot(vp_hyde1, aes(x=L1, y=value,fill=stat))+geom_violin()+labs(x="Triplets", y = "Value")+scale_fill_manual(values=c("Grey", "goldenrod2"),name="",labels=c("D",expression(gamma)))+ geom_boxplot(width=0.05,position=position_dodge(0.9),outlier.size=-1)+stat_summary(fun.y=median, geom="point", size=2, color="red",position=position_dodge(0.9))




#Total tree test
total=read.table("hyde_all_tri-out_noephemera.txt",header=T)
total=total[complete.cases(total), ]
total$D=(total$ABBA-total$ABAB)/(total$ABBA+total$ABAB)
total=total[total$Gamma > 0 & total$Gamma < 1, ]
total$Pvalue=p.adjust(total$Pvalue,method="bonferroni")

dat=total[,c("Gamma","D","Pvalue")]
g1=ggplot(dat)+geom_histogram(aes(Gamma),bins=500,col="black")+labs(x = expression(gamma),y="N quartets")
g2=ggplot(dat)+geom_histogram(aes(D),bins=500,col="black")+labs(x = "D",y="N quartets")
dat$density = get_density(dat$Gamma, dat$D, n = 100)
g3=ggplot(dat) + geom_point(aes(Gamma, D, color = density),size=0.1)+scale_color_gradientn(colors = cividis(30))+labs(x = expression(gamma),y="D")
g4=ggplot(dat) + geom_point(aes(Gamma, D, color = Pvalue),size=0.1)+scale_color_gradientn(colors = cividis(30))+labs(x = expression(gamma),y="D")

dat=dat[dat$Pvalue<0.05,c("Gamma","D","Pvalue")]
g5=ggplot(dat)+geom_histogram(aes(Gamma),bins=500,col="black")+labs(x = expression(gamma),y="N quartets")
g6=ggplot(dat)+geom_histogram(aes(D),bins=500,col="black")+labs(x = "D",y="N quartets")
dat$density = get_density(dat$Gamma, dat$D, n = 100)
g7=ggplot(dat) + geom_point(aes(Gamma, D, color = density),size=0.1)+scale_color_gradientn(colors = cividis(30))+labs(x = expression(gamma),y="D")
g8=ggplot(dat) + geom_point(aes(Gamma, D, color = Pvalue),size=0.1)+scale_color_gradientn(colors = cividis(30))+labs(x = expression(gamma),y="D")

size 14.348624  4.834862
grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=4,nrow=2)
quartz.save("hyde_total", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)

#PHATE
total_mult=total[total$Pvalue<0.05,]
total_phate=total_mult[,c(7:21)]/apply(total_mult[,c(7:21)],1,sum)
phate_out=phate(total_phate)
dat1=data.frame(phate_out$embedding)
dat1$Gamma=total_mult$Gamma
dat1$D=total_mult$D

g9=ggplot(dat1) + geom_point(aes(PHATE1, PHATE2, color = Gamma),size=0.01)+scale_color_gradientn(colors = cividis(30),name=expression(gamma))+theme_classic()
g10=ggplot(dat1) + geom_point(aes(PHATE1, PHATE2, color = D),size=0.01)+scale_color_gradientn(colors = cividis(30))+theme_classic()
grid.arrange(g9,g10,ncol=2)
quartz.save("phate_total", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)