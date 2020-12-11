library('ape')
library('phytools')
library('phangorn')
library('Rtsne')
library('MASS')
library('gridExtra')
library('ggplot2')
library('pals')
library('reshape2')
library("svMisc")
#library("VennDiagram")
library("MCMCtreeR")
library("repr")
setwd("~/Desktop/data/odonata_introg_data/")

################################################################ HyDe ######################################################################## 
get_density = function(x, y, ...) 
{
  dens = MASS::kde2d(x, y, ...)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}



###Wilcoxon Rank sum test for Ds and Gammas 
get_wilx_hyde=function(m,stats)
{
    order_mean=mean(m[,stats])
    m_out=c()
    for(i in c("Suborder","Superfamily_between","Superfamily_within"))
    {
       
        for (cl in unique(m[,i]))
        {
            
            my_mean=mean(m[m[,i]==cl,stats])
            my_pt=wilcox.test(m[m[,i]==cl,stats],m[m[,i]!=cl,stats],alternative="two.sided")$p.value
            m_out=rbind(m_out,c(i,cl,my_mean,my_pt,my_pt<0.05,order_mean<my_mean))
        }
            
    }
    m_out=as.data.frame(m_out)
    names(m_out)=c("rank","taxon","average","p_two","sig","greater")
    m_out=m_out[m_out$taxon!="RANDOM",]
    return(m_out)
        
}
###Fisher exact test for excess of significant triplets 
get_fet_hyde=function(m)
{
    all_i=table(m$introgressionid)
    m_out=c()
    for(i in c("Suborder","Superfamily_between","Superfamily_within"))
    {
       
        for (cl in unique(m[,i]))
        {
            
            cl_i=table(m[m[,i]==cl,"introgressionid"])
            m_fet=cbind(all_i,cl_i)
            pval=fisher.test(m_fet)$p.value
            
            m_out=rbind(m_out,c(i,cl,all_i[1],all_i[2],cl_i[1],cl_i[2],pval,pval<0.05,cl_i[1]/cl_i[2],(cl_i[1]/cl_i[2])>(all_i[1]/all_i[2])))
        }
            
    }
    m_out=as.data.frame(m_out)
    names(m_out)=c("rank","taxon","Total_itrog_count","Total_noitrog_count","Clade_itrog_count","Clade_noitrog_count","P_FET","Significant(<0.05)","Ratio","Excess")
    m_out=m_out[m_out$taxon!="RANDOM",]
    return(m_out) 
    
}    





###Truplet consistency with species topology/ D statistic/ Gamma/ putative introgressing taxa
good_tri=function(table_row,phy)
{
    P1=unlist(table_row$P1)
    P2=unlist(table_row$P2)
    H=unlist(table_row$Hybrid)
    tri=keep.tip(phy,c(P1,H,P2))
    P1H=any(c(is.monophyletic(tri,c(P1,H)),is.monophyletic(tri,c(H,P1))))   
    P2H=any(c(is.monophyletic(tri,c(P2,H)),is.monophyletic(tri,c(H,P2))))
    if (P1H)
    {
        trip_good=TRUE
        intr=sort(c(P2,H))
    }else if (P2H)
    {
        trip_good=TRUE
        intr=sort(c(P1,H))
    }else{
      
        trip_good=FALSE
        intr=c(NaN,NaN)
    }
    return(data.frame(trip_consistent=trip_good,i1=intr[1],i2=intr[2]))       
}    

tt=as.character(read.table("taxa_map.txt")$V1)
phy=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")
phy_mcmc=readMCMCtree("FigTree4long.tre")
phy_dated=phy_mcmc$apePhy
Zygoptera=tt[1:48]
Anisoptera=tt[50:83]
Anisozygoptera=tt[49]
Epiprocta=tt[49:83]
Lestoidea=tt[44:48]
Calopterygoidea=tt[1:19]
Coenagrionoidea=tt[20:42]
Aeshnoidea=tt[50:63]
Cordulegastroidea=tt[64:69]
Libelluloidea=tt[70:83]
Platystictidae=tt[43]
Aeshnidae=tt[50:56]
Gomphidae_Petaluridae=tt[57:63]
Libellulidae=tt[75:83]


total=read.table("hyde_all_tri-out_noephemera.txt",header=T,stringsAsFactors=F)
total=total[complete.cases(total$Gamma),]
total=total[total$Gamma<=1 & total$Gamma>=0, ]
names(total)[6]="Gamma01"
total$Gamma=ifelse(total$Gamma01>0.5,1-total$Gamma01,total$Gamma01)
total$polarized=ifelse(total$ABBA>total$AABB & total$AABB<total$ABAB, "FALSE",ifelse(total$ABBA<total$AABB & total$ABAB>total$ABBA, "FALSE","TRUE"))
total=total[total$polarized==TRUE,]
total$P1=unlist(lapply(lapply(strsplit(as.vector(total$P1),"_"),"[",c(2,3)),paste,collapse="_"))
total$Hybrid=unlist(lapply(lapply(strsplit(as.vector(total$Hybrid),"_"),"[",c(2,3)),paste,collapse="_"))
total$P2=unlist(lapply(lapply(strsplit(as.vector(total$P2),"_"),"[",c(2,3)),paste,collapse="_"))

#Triplet age
tree_h=nodeheight(phy_dated,node=85)
total$t_P1H=tree_h-apply(cbind(unlist(total[,1]),unlist(total[,2])),1,findMRCA,tree=phy_dated,type="height")
total$t_P2H=tree_h-apply(cbind(unlist(total[,2]),unlist(total[,3])),1,findMRCA,tree=phy_dated,type="height")
total$t_P1P2=tree_h-apply(cbind(unlist(total[,1]),unlist(total[,3])),1,findMRCA,tree=phy_dated,type="height")

#D / Chisq test for D significance
total$D1=(total$ABBA-total$ABAB)/(total$ABBA+total$ABAB)
total$D2=(total$ABAB-total$AABB)/(total$ABAB+total$AABB)
total$chiPd1=unlist(lapply(apply(total[,c("ABBA","ABAB")],1,chisq.test),"[[","p.value"))
total$chiPd2=unlist(lapply(apply(total[,c("ABAB","AABB")],1,chisq.test),"[[","p.value"))
total$D=abs(ifelse(total$ABBA>total$AABB,total$D2,total$D1))
total$chiPd=ifelse(total$ABBA>total$AABB,total$chiPd2,total$chiPd1)
total$PvalueD=p.adjust(total$chiPd,method="bonferroni")
total$Pvalue=p.adjust(total$Pvalue,method="bonferroni")


s=c()
for (i in 1:nrow(total))
{
    progress(i,nrow(total))
    v=good_tri(total[i,],phy)
    s=rbind(s,v)
}
total=cbind(total,s)
total=total[total$trip_consistent==TRUE,]


#Taxon assignment 
total$P1subo=ifelse(total$P1 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P1 %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$Hybridsubo=ifelse(total$Hybrid %in% Zygoptera,"Zygoptera",
                    ifelse(total$Hybrid %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$P2subo=ifelse(total$P2 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P2 %in% Anisoptera,"Anisoptera","Anisozygoptera"))
    

total$P1sub=ifelse(total$P1 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P1 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P1 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P1 %in% Platystictidae,"Platystictoidea",
                   ifelse(total$P1 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P1 %in% Platystictidae,"Platystictoidea",       
                   ifelse(total$P1 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P1 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P1 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P1 %in% Anisozygoptera,"Epiophlebiidae","RANDOM"))))))))))

total$Hybridsub=ifelse(total$Hybrid %in% Lestoidea,"Lestoidea",
                   ifelse(total$Hybrid %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$Hybrid %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$Hybrid %in% Platystictidae,"Platystictoidea",
                   ifelse(total$Hybrid %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$Hybrid %in% Platystictidae,"Platystictoidea",       
                   ifelse(total$Hybrid %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$Hybrid %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$Hybrid %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$Hybrid %in% Anisozygoptera,"Epiophlebiidae","RANDOM"))))))))))


total$P2sub=ifelse(total$P2 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P2 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P2 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P2 %in% Platystictidae,"Platystictoidea",
                   ifelse(total$P2 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P1 %in% Platystictidae,"Platystictoidea",
                   ifelse(total$P2 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P2 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P2 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P2 %in% Anisozygoptera,"Epiophlebiidae","RANDOM"))))))))))

#Ephiophlebia only between Zygoptera and Anisoptera 
#nodups=!apply(t(apply(total[,c("P1sub","Hybridsub","P2sub")],1,duplicated)),1,any)

total$Order="Odonata"

total$Suborder=ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Zygoptera),1,all) ,"Zygoptera",
                      ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Epiprocta),1,all),"Epiprocta","RANDOM"))
                     
total$Superfamily_between=ifelse(total$Hybrid %in% Anisozygoptera,"Anisozygoptera",
                         ifelse(total$Hybrid %in% Lestoidea & apply(apply(total[,c("P1","P2")],2,"%in%",c(Epiprocta,Calopterygoidea,Coenagrionoidea,Platystictidae)),1,all),"Lestoidea",
                         ifelse(total$Hybrid %in% Platystictidae & apply(apply(total[,c("P1","P2")],2,"%in%",c(Lestoidea,Calopterygoidea,Coenagrionoidea)),1,all),"Platystictoidea",
                         ifelse(total$Hybrid %in% Coenagrionoidea & apply(apply(total[,c("P1","P2")],2,"%in%",c(Platystictidae,Calopterygoidea,Lestoidea)),1,all),"Coenagrionoidea",
                         ifelse(total$Hybrid %in% Calopterygoidea & apply(apply(total[,c("P1","P2")],2,"%in%",c(Platystictidae,Coenagrionoidea,Lestoidea)),1,all),"Calopterygoidea",
                         ifelse(total$Hybrid %in% Aeshnidae & apply(apply(total[,c("P1","P2")],2,"%in%",c(Anisozygoptera,Gomphidae_Petaluridae,Cordulegastroidea,Libelluloidea)),1,all),"Aeshnoidea (Aeshnidae)",
                         ifelse(total$Hybrid %in% Gomphidae_Petaluridae & apply(apply(total[,c("P1","P2")],2,"%in%",c(Anisozygoptera,Aeshnidae,Cordulegastroidea,Libelluloidea)),1,all),"Aeshnoidea (Gomphidae+Petaluridae)",
                         ifelse(total$Hybrid %in% Cordulegastroidea & apply(apply(total[,c("P1","P2")],2,"%in%",c(Anisozygoptera,Aeshnoidea,Libelluloidea)),1,all),"Cordulegastroidea",
                         ifelse(total$Hybrid %in% Libelluloidea & apply(apply(total[,c("P1","P2")],2,"%in%",c(Anisozygoptera,Aeshnoidea,Cordulegastroidea)),1,all),"Libelluloidea","RANDOM")))))))))       
                                
                                                                
total$Superfamily_within=ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Aeshnidae),1,all),"Aeshnoidea (Aeshnidae)",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Gomphidae_Petaluridae),1,all),"Aeshnoidea (Gomphidae+Petaluridae)",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(total[,c("P1","P2","Hybrid")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM")))))))       
                              

total$introgressionid=ifelse(total$PvalueD<0.05 & total$Pvalue<10^-6,"Introgression","None") 
#Calculate density for scatter plot



write.csv(total,"total_hyde.csv",quote=F,row.names=F)
total=read.csv("total_hyde.csv",stringsAsFactors=F)



#Color scheme
sunset=c("#364B9A" ,"#4A7BB7", "#6EA6CD", "#98CAE1" ,"#C2E4EF" ,"#EAECCC", "#FEDA8B" ,"#FDB366", "#F67E4B", "#DD3D2D", "#A50026")

#D and Gamma distributions Violin plots
total_ord=melt(total[,c("D","Gamma","Pvalue","Order","Suborder","Superfamily_between","Superfamily_within","PvalueD")],value.name="taxon",id=c("D","Gamma","Pvalue","PvalueD"))
total_ord=total_ord[total_ord$taxon!="RANDOM",]
total_ord$variable=replace(as.character(total_ord$variable),as.character(total_ord$variable)=="Superfamily_within","Superfamily within")
total_ord$variable=replace(as.character(total_ord$variable),as.character(total_ord$variable)=="Superfamily_between","Superfamily between")
total_ord$variable_f=factor(total_ord$variable, levels=c("Order","Suborder","Superfamily between","Superfamily within"))
total_ord$taxon=factor(total_ord$taxon, levels=levels(factor(total_ord$taxon))[c(10,7,12,3,1,2,6,9,8,11,4,5)])


###############Plotting D
a1=ggplot(total_ord[total_ord$PvalueD<0.05,], aes(x=taxon, y=D))+geom_violin(lwd=0.1,fill='olivedrab3')+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+theme(axis.text.x = element_text(size = 7,angle=15,hjust = 1),plot.margin=unit(c(0,0.1,0,0.1), "cm"))+ylab("D")+xlab("")+ggtitle("A")

###############Plotting Gamma
a2=ggplot(total_ord[total_ord$Pvalue<10^-6,], aes(x=taxon, y=Gamma))+geom_violin(lwd=0.1,fill='cornflowerblue')+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun=mean, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+theme(axis.text.x = element_text(size = 7,angle=15,hjust = 1),plot.margin=unit(c(0,0.1,0,0.1), "cm"))+ylab(expression(gamma))+xlab("")+ggtitle("B")

###############Plotting D vs Gamma
total_s=total[total$introgressionid=="Introgression",c("Gamma","D")]
total_s$density = get_density(total_s$Gamma, total_s$D, n = 300)

a4=ggplot(total_s,aes(Gamma, D, color = density)) + geom_point(size=0.4,stroke=0)+geom_smooth(method = "auto", size = 0.5,color="black")+scale_color_gradientn(colors = sunset,guide=guide_colorbar(barwidth = 6,barheight =0.5,label.theme = element_text(size=8)),name="")+labs(x = expression(gamma),y="D")+theme(legend.position=c(0.5,1.05) ,legend.direction="horizontal",legend.background = element_blank())+ggtitle("D")



###############Plotting Counts
total_p=melt(total[,c("introgressionid","Order","Suborder","Superfamily_between","Superfamily_within")],id="introgressionid",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily_between","Superfamily_within"))
total_p$variable_f=replace(as.character(total_p$variable_f),as.character(total_p$variable_f)=="Superfamily_within","Superfamily within")
total_p$variable_f=replace(as.character(total_p$variable_f),as.character(total_p$variable_f)=="Superfamily_between","Superfamily between")
total_p$taxon=factor(total_p$taxon, levels=levels(factor(total_p$taxon))[c(10,7,12,3,1,2,6,9,8,11,4,5)])


a3=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=introgressionid))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5),size=2.5)+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1),legend.position=c(0.5,1.4) ,legend.direction="horizontal",legend.background = element_blank())+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("wheat","grey50"),name="")+ggtitle("C")


#tSNE
totalsig=total[total$introgressionid=="Introgression" & total$D <0.25,c("Gamma","D")]
total_dr=total[total$PvalueD<0.05 & total$D <0.25 & total$Pvalue<10^-6,c(7:21)]/apply(total[total$PvalueD<0.05 & total$D<0.25 & total$Pvalue<10^-6,c(7:21)],1,sum)
tsne_out=Rtsne(total_dr)
dat1=data.frame(tsne_out$Y)
names(dat1)=c("tSNE1","tSNE2")
dat1=cbind(dat1,totalsig)

###############Plotting tSNE D
a5=ggplot(dat1) + geom_point(aes(tSNE1, tSNE2, color = D),size=0.35,stroke=0)+scale_color_gradientn(colors = sunset,name="D",guide=guide_colorbar(barwidth = 6,barheight =0.5,label.theme = element_text(size=8)))+theme_classic()+theme(legend.position=c(0.5,1.05) ,legend.direction="horizontal",strip.background = element_rect(colour = "white", fill = "white"),legend.background = element_blank())+ggtitle("E")

###############Plotting tSNE D Gamma
a6=ggplot(dat1) + geom_point(aes(tSNE1, tSNE2, color = Gamma),size=0.35,stroke=0)+scale_color_gradientn(colors = sunset,name=expression(gamma),guide=guide_colorbar(barwidth = 6,barheight =0.5,label.theme = element_text(size=8)))+theme_classic()+theme(legend.position=c(0.5,1.05) ,legend.direction="horizontal",legend.background = element_blank())+ggtitle("F")

lay = rbind(c(1,1,1),
             c(2,2,2),
           c(3,3,3),
            c(4,5,6))


f=arrangeGrob(a1,a2,a3,a4,a5,a6,layout_matrix = lay)
ggsave("Fig3.png",f, width = 7, height = 9.5)
ggsave("Fig3.pdf",f, width = 7, height = 9.5)


#quartz.save("Fig3.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
#quartz.save("Fig3.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)

#####################HyDe tests#####################
#D for old vs young divergencies
old=total[total$PvalueD<0.05 & total$Superfamily_between!="RANDOM","D"]
new=total[total$PvalueD<0.05 & total$Superfamily_within!="RANDOM","D"]
mean(old)
mean(new)
wilcox.test(old,new)
very_new=total[total$PvalueD<0.05 & (total$Superfamily_within %in% c("Calopterygoidea","Coenagrionoidea","Libelluloidea") ),"D"]
mean(very_new)
wilcox.test(old,very_new)

#D odonata vs. all
get_wilx_hyde(total[total$PvalueD<0.05,],"D")

#Gamma odonata vs. all
get_wilx_hyde(total[total$Pvalue<10^-6,],"Gamma")

#Excess of introgression triplets odonata vs. all
get_fet_hyde(total)


#Correlation test between Gamma and D
a=total[total$introgressionid=="Introgression",c("Gamma","D")
cor.test(a[,1],a[,2])







#Introgressing pair

introg=data.frame(t(apply(cbind(unlist(lapply(lapply(strsplit(as.vector(total$P1),"_"),"[",c(2,3)),paste,collapse="_")),unlist(lapply(lapply(strsplit(as.vector(total$Hybrid),"_"),"[",c(2,3)),paste,collapse="_"))),1,sort)))
names(introg)=c("i1","i2")
i1i2=as.vector(apply(introg,1,paste,collapse="_"))


total_hyde=cbind(total,introg)
total_hyde$i1i2=i1i2

################################################################ Dfoil ########################################################################

get_wilx_t_dfoil=function(m,stats)
{
    m_out=c()
    for(i in c("Suborder","Superfamily_within","Superfamily_between"))
    {
       
        for (cl in unique(m[,i]))
        {
            
            my_mean=mean(unlist(m[m[,i]==cl,stats]))
            my_pt=wilcox.test(unlist(m[m[,i]==cl,stats]),unlist(m[,stats]),alternative="two.sided")$p.value
            my_t=t.test(m[m[,i]==cl,stats],mu=0)$p.value
            m_out=rbind(m_out,c(i,cl,my_mean,my_pt,my_t,my_pt<0.05,my_t<0.05))
            
        }
            
    }
    m_out=as.data.frame(m_out)
    names(m_out)=c("Rank","Taxon","Average_Dfoil","P_wilx","P_t","Significant_wilx","Significant_t")
    m_out=m_out[m_out$Taxon!="RANDOM",]
    return(m_out)
        
}

get_fet_dfoil=function(m)
{
    all_i=table(m$introgressionid)
    noi=all_i[names(all_i)=="None"]
    intr=sum(all_i[names(all_i)!="None"])
    all_i=c(intr,noi)
    
    m_out=c()
    for(i in c("Suborder","Superfamily_between","Superfamily_within"))
    {
       
        for (cl in unique(m[,i]))
        {
            
            cl_i=table(m[m[,i]==cl,"introgressionid"])
           
            if (length(cl_i) > 1)
            {
                noir=cl_i[names(cl_i)=="None"]
                intr=sum(cl_i[names(cl_i)!="None"])
                cl_i=c(intr,noir)
                m_fet=cbind(all_i,cl_i)
                pval=fisher.test(m_fet)$p.value
                m_out=rbind(m_out,c(i,cl,all_i[1],all_i[2],cl_i[1],cl_i[2],pval,pval<0.05,(all_i[1]/all_i[2])<(cl_i[1]/cl_i[2])))
            }    
            
        }
            
    }
    m_out=as.data.frame(m_out)
    names(m_out)=c("rank","taxon","Total_itrog_count","Total_noitrog_count","Clade_itrog_count","Clade_noitrog_count","P_FET","Significant(<0.05)","Excess")
    m_out=m_out[m_out$taxon!="RANDOM",]
    return(m_out) 
    
}    

        
        
        
        
        
names_v=c("P1","P2","P3","P4","Out","chrom1","position", "AAAAA" , "AAABA" , "AABAA" , "AABBA" , "ABAAA" , "ABABA" , "ABBAA" , "ABBBA" , "BAAAA" , "BAABA" , "BABAA" , "BABBA" , "BBAAA" ,"BBABA","BBBAA","BBBBA",'chromdup','coord','total','dtotal','T12','T34','T1234','DFO_left','DFO_right','DFO_total','DFO_stat','DFO_chisq','DFO_Pvalue','DIL_left','DIL_right','DIL_total','DIL_stat','DIL_chisq','DIL_Pvalue','DFI_left','DFI_right','DFI_total','DFI_stat','DFI_chisq','DFI_Pvalue','DOL_left','DOL_right','DOL_total','DOL_stat','DOL_chisq','DOL_Pvalue','introgression','introgna','intrognone','introg13','introg14','introg23','introg24','introg31','introg41','introg32','introg42','introg123','introg124')

tt=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")
Anisozygoptera=c("Epiophlebia_superstes")
Zygoptera=c(extract.clade(tt,88)$tip.label)
Anisoptera=c(extract.clade(tt,136)$tip.label)
Lestoidea=c(extract.clade(tt,131)$tip.label)
Calopterygoidea=extract.clade(tt,91)$tip.label
Coenagrionoidea=extract.clade(tt,109)$tip.label
Cordulegastroidea=extract.clade(tt,151)$tip.label
Libelluloidea=extract.clade(tt,156)$tip.label
Platystictidae="Protosticta_beaumonti"
Aeshnidae=extract.clade(tt,137)$tip.label
Gomphidae_Petaluridae=extract.clade(tt,144)$tip.label
Epiprocta=c(Anisoptera,Anisozygoptera)
phy_mcmc=readMCMCtree("FigTree4long.tre")
phy_dated=phy_mcmc$apePhy



total=read.csv("dfoil_results.txt",stringsAsFactors=FALSE)
names(total)=names_v

#Triplet age
tree_h=nodeheight(phy_dated,node=85)
total$t_P1P2=tree_h-apply(cbind(unlist(total[,1]),unlist(total[,2])),1,findMRCA,tree=phy_dated,type="height")
total$t_P3P4=tree_h-apply(cbind(unlist(total[,3]),unlist(total[,4])),1,findMRCA,tree=phy_dated,type="height")        
        
#T12<T34

total=total[total$t_P3P4>total$t_P1P2,]
        


total$P1subo=ifelse(total$P1 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P1 %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$P2subo=ifelse(total$P2 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P2 %in% Anisoptera,"Anisoptera","Anisozygoptera"))


total$P3subo=ifelse(total$P3 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P3 %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$P4subo=ifelse(total$P4 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P4 %in% Anisoptera,"Anisoptera","Anisozygoptera"))


total$P1sub=ifelse(total$P1 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P1 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P1 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P1 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P1 %in% Aeshnidae,"Aeshnoidea(Aeshnidae)",
                   ifelse(total$P1 %in% Gomphidae_Petaluridae,"Aeshnoidea(Gomphidae+Petaluridae)",
                   ifelse(total$P1 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P1 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P1 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))

total$P3sub=ifelse(total$P3 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P3 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P3 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P3 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P3 %in% Aeshnidae,"Aeshnoidea(Aeshnidae)",
                   ifelse(total$P3 %in% Gomphidae_Petaluridae,"Aeshnoidea(Gomphidae+Petaluridae)",
                   ifelse(total$P3 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P3 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P3 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))


total$P2sub=ifelse(total$P2 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P2 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P2 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P2 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P2 %in% Aeshnidae,"Aeshnoidea(Aeshnidae)",
                   ifelse(total$P2 %in% Gomphidae_Petaluridae,"Aeshnoidea(Gomphidae+Petaluridae)",
                   ifelse(total$P2 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P2 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P2 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))


total$P4sub=ifelse(total$P4 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P4 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P4 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P4 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P4 %in% Aeshnidae,"Aeshnoidea(Aeshnidae)",
                   ifelse(total$P4 %in% Gomphidae_Petaluridae,"Aeshnoidea(Gomphidae+Petaluridae)",
                   ifelse(total$P4 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P4 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P4 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))

total$Order="Odonata"

total$Suborder=ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Zygoptera),1,all) ,"Zygoptera",
                      ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Anisoptera),1,all) ,"Anisoptera","RANDOM"))
                     

total$Superfamily_between=ifelse(apply(apply(total[,1:4],2,"%in%",c("Epiophlebia_superstes")),1,any),"1",
                         ifelse(apply(apply(total[,1:2],2,"%in%",Lestoidea),1,all) | apply(apply(total[,3:4],2,"%in%",Lestoidea),1,all),"5",
                         ifelse((apply(apply(total[,1:2],2,"%in%",Coenagrionoidea),1,all) & apply(apply(total[,3:4],2,"%in%",Calopterygoidea),1,all)) | (apply(apply(total[,1:2],2,"%in%",Calopterygoidea),1,all) & apply(apply(total[,3:4],2,"%in%",Coenagrionoidea),1,all)),"6",
                         ifelse((apply(apply(total[,1:2],2,"%in%",Aeshnidae),1,all) & apply(apply(total[,3:4],2,"%in%",c(Gomphidae_Petaluridae,Cordulegastroidea,Libelluloidea)),1,all)) | (apply(apply(total[,1:2],2,"%in%",c(Gomphidae_Petaluridae,Cordulegastroidea,Libelluloidea)),1,all) & apply(apply(total[,3:4],2,"%in%",Aeshnidae),1,all)),"2",
                         ifelse((apply(apply(total[,1:2],2,"%in%",Gomphidae_Petaluridae),1,all) & apply(apply(total[,3:4],2,"%in%",c(Cordulegastroidea,Libelluloidea)),1,all)) | (apply(apply(total[,1:2],2,"%in%",c(Cordulegastroidea,Libelluloidea)),1,all) & apply(apply(total[,3:4],2,"%in%",Gomphidae_Petaluridae),1,all)),"3",
                         ifelse((apply(apply(total[,1:2],2,"%in%",Cordulegastroidea),1,all) & apply(apply(total[,3:4],2,"%in%",Libelluloidea),1,all)) | (apply(apply(total[,1:2],2,"%in%",Libelluloidea),1,all) & apply(apply(total[,3:4],2,"%in%",Cordulegastroidea),1,all)),"4","RANDOM"))))))     
                         
            
        

total$Superfamily_within=ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Aeshnidae),1,all),"Aeshnoidea(Aeshnidae)",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Gomphidae_Petaluridae),1,all),"Aeshnoidea(Gomphidae+Petaluridae)",       
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM")))))))


#Subsampling Anisozygoptera
total_noaniszyg= total[ !apply(apply(total[,1:4],2,"%in%",c("Epiophlebia_superstes")),1,any),]
total_aniszyg= total[ apply(apply(total[,1:4],2,"%in%",c("Epiophlebia_superstes")),1,any),][sample(1:38352,size=4000),] 
total_sub=rbind(total_noaniszyg,total_aniszyg)   
total_sub$DFO_Pvalue=p.adjust(total_sub$DFO_Pvalue,method="fdr")
total_sub$DIL_Pvalue=p.adjust(total_sub$DIL_Pvalue,method="fdr")
total_sub$DFI_Pvalue=p.adjust(total_sub$DFI_Pvalue,method="fdr")
total_sub$DOL_Pvalue=p.adjust(total_sub$DOL_Pvalue,method="fdr")
P_adj=apply(apply(total_sub[,c("DFO_Pvalue","DIL_Pvalue","DFI_Pvalue","DOL_Pvalue")],2,"<",0.05),1,all)
total_sub$introgression_corr=ifelse(P_adj==FALSE,"none",total_sub$introgression)         
total_sub$introgressionid=ifelse(total_sub$introgression_corr=="none","None",
                             ifelse(total_sub$introgression_corr=="123" | total_sub$introgression_corr=="124","Ancestral","Inter-group"))        
        
        
#Dfoil violin plots        
total_ord=melt(total_sub[total_sub$introgressionid!="None" ,c("DFO_stat","DIL_stat","DFI_stat","DOL_stat","Order","Superfamily_within","Suborder","Superfamily_between")],measure.vars = c("DFO_stat","DIL_stat","DFI_stat","DOL_stat"),value.name="stat")      
total_ord=melt(total_ord,id=c("variable","stat"),value.name="taxon",variable.name="class")
total_ord$class=replace(as.character(total_ord$class),as.character(total_ord$class)=="Superfamily_between","Superfamily between")
total_ord$class=replace(as.character(total_ord$class),as.character(total_ord$class)=="Superfamily_within","Superfamily within")
total_ord$variable_f=factor(total_ord$class, levels=c("Order","Suborder","Superfamily between","Superfamily within"))
total_ord$variable=ifelse(total_ord$variable=="DFO_stat","D[FO]",ifelse(total_ord$variable=="DIL_stat","D[IL]",ifelse(total_ord$variable=="DFI_stat","D[FI]","D[OL]")))


###############Plotting 
g1=ggplot(total_ord[total_ord$taxon!="RANDOM" ,], aes(x=taxon, y=stat))+geom_violin(lwd=0.1,fill="plum2")+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun=median, geom="point", size=2,position=position_dodge(0.9))+geom_boxplot(width=0.01,outlier.size=-1,position=position_dodge(0.9))+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1))+ylab(expression(D[FOIL]))+xlab("")+theme(legend.position=c(0.5,1.4),legend.direction="horizontal",legend.background = element_blank())+ggtitle("B") 




#Proportions
total_p=melt(total_sub[,c("introgressionid","Order","Superfamily_within","Suborder","Superfamily_between")],id="introgressionid",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="Superfamily_between","Superfamily between")
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="Superfamily_within","Superfamily within")        
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily between","Superfamily within"))

###############Plotting
g2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=introgressionid))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5),size=2)+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1),legend.position=c(0.5,1.4),legend.direction="horizontal",legend.background = element_blank())+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("rosybrown2","wheat","grey50"),name="")+ggtitle("C")

lay = rbind(c(1,1,1),
             c(2,2,2))
           


f=arrangeGrob(g1,g2,layout_matrix = lay)
ggsave("Fig4.png",f, width = 7, height = 6)


        
#####################DFoil tests#####################        
        
#Dfoil odonata vs. all        
get_wilx_t_dfoil(total_sub[total_sub$introgressionid!="None",],c("DFO_stat","DIL_stat","DFI_stat","DOL_stat"))

#Fet odonata vs. all         
get_fet_dfoil(total_sub)
    
#Introgressing pair
#Identify introgressing taxa pair from Dfoil result
get_intropair_dfoil=function(m)
{
    pair_v=rbind()
    for (i in 1:nrow(m))
    {
        progress(i,nrow(m))
        if(m[i,"introgression"]!="none" & m[i,"introgression"]!="123" & m[i,"introgression"]!="124")
        {
            pair=sort(as.character(unlist(m[i,c("P1","P2","P3","P4")][sort(as.numeric(unlist(strsplit(as.character(m[i,"introgression"]),split=""))))])))
            pair_v=rbind(pair_v,pair)
        } else if (m[i,"introgression"]=="none") {
            pair=sort(as.character(unlist(m[i,c(sample(c("P1","P2"),1),sample(c("P3","P4"),1))])))
            pair_v=rbind(pair_v,pair)
        } else {    
            pair=sort(as.character(unlist(m[i,c("P1","P2","P3","P4")][sort(as.numeric(unlist(strsplit(as.character(m[i,"introgression"]),split=""))))[c(sample(1:2,1),3)]])))
            pair_v=rbind(pair_v,pair)
        }   
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1","i2")
    rownames(pair_v) = c()
    return(pair_v)
}

intropairs=get_intropair_dfoil(total)
i1i2=as.vector(apply(intropairs,1,paste,collapse="_"))
total_dfoil=cbind(total,intropairs)
total_dfoil$i1i2=i1i2


        
        
        
        
        
        
################################################################ Chi-square/BLT ############################################################################        
tt=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")
Anisozygoptera=c("Epiophlebia_superstes")
Zygoptera=c(extract.clade(tt,88)$tip.label)
Anisoptera=c(extract.clade(tt,136)$tip.label)
Lestoidea=c(extract.clade(tt,131)$tip.label)
Calopterygoidea=extract.clade(tt,91)$tip.label
Coenagrionoidea=extract.clade(tt,109)$tip.label
Cordulegastroidea=extract.clade(tt,151)$tip.label
Libelluloidea=extract.clade(tt,156)$tip.label
Platystictidae="Protosticta_beaumonti"
Aeshnidae=extract.clade(tt,137)$tip.label
Gomphidae_Petaluridae=extract.clade(tt,144)$tip.label
Epiprocta=c(Anisoptera,Anisozygoptera)
c234=c(Aeshnidae,Gomphidae_Petaluridae,Cordulegastroidea,Libelluloidea)
c34=c(Gomphidae_Petaluridae,Cordulegastroidea,Libelluloidea)
c4=c(Cordulegastroidea,Libelluloidea)       
c86=c(Platystictidae,Calopterygoidea,Coenagrionoidea)        
c6=c(Calopterygoidea,Coenagrionoidea)        
        
total=read.csv("chi_sq_all.txt",stringsAsFactors=FALSE,header=F)        
names_vb=c("clade","P1out","P2out","P3out","CountP1","CountP2","CountP3","congruent","PvalueChi","meanT_concord","meanT_discord1","meanT_discord2","PvalueWCOMC1","PvalueWCOMC2","PvalueWC1C2")
names(total)=names_vb       
#Chi square
total$PvalueChi=p.adjust(total$PvalueChi,method="fdr")        
total$PvalueChiExtr=p.adjust(unlist(lapply(apply(total[,c("CountP1","CountP2","CountP3")],1,chisq.test),"[[","p.value")),method="fdr")        
total$PvalueWC1C2=p.adjust(total$PvalueWC1C2,method="fdr")
        
total$introgressionid=ifelse(total$PvalueChiExtr > 0.05,"Extreme ILS",ifelse(total$PvalueChi>0.05 | total$PvalueWC1C2 >0.05 | total$PvalueWCOMC2  >0.05 | total$PvalueWCOMC1  >0.05 ,"ILS",ifelse(total$PvalueChi<0.05 & total$PvalueWC1C2 < 0.05 & (total$meanT_concord < total$meanT_discord2)  & (total$meanT_discord2 < total$meanT_discord1) ,"Introgression+ILS","ILS")))        
        
        
total$Order="Odonata"

total$Suborder=ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Zygoptera),1,all) ,"Zygoptera",
                      ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Epiprocta),1,all) ,"Epiprocta","RANDOM"))
                     

total$Superfamily_between=ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Anisozygoptera),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Anisoptera),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Zygoptera),1,any) ,"1",
                           ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Aeshnidae),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",c34),1,any),"2",
                           ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Gomphidae_Petaluridae),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",c4),1,any),"3",       
                           ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Cordulegastroidea),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Libelluloidea),1,any),"4",      
                           ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Lestoidea),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",c86),1,any),"5",
                           ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Platystictidae),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",c6),1,any),"8",
                           ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Coenagrionoidea),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Calopterygoidea),1,any),"6",
                           ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Anisozygoptera),1,any) & apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Anisoptera),1,any),"7","RANDOM"))))))))      
        

total$Superfamily_within=ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Aeshnidae),1,all),"Aeshnoidea(Aeshnidae)",
                       ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Gomphidae_Petaluridae),1,all),"Aeshnoidea(Gomphidae+Petaluridae)",       
                       ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(total[,c("P1out","P2out","P3out")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM")))))))        
total=total[total$congruent==T,]        
        
        
        
#Proportions
total_p=melt(total[,c("introgressionid","Order","Superfamily_within","Suborder","Superfamily_between")],id="introgressionid",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="Superfamily_between","Superfamily between")
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="Superfamily_within","Superfamily within")        
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily between","Superfamily within"))

###############Plotting
ch1=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=introgressionid))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5),size=2)+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1),legend.position=c(0.5,1.25),legend.direction="horizontal",legend.background = element_blank())+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("gold","grey50","wheat","red"),name="")+ggtitle("B")        
        

################Branch Lengths         
bl=read.table("all_brls.csv",stringsAsFactors=FALSE,header=F)          
names_bl=c("clade_name","P1","P2","P3","brl1","brl2","trl","brl_out","brl_int","root_tip","topo","topo_con_sp","good_trip","tripl","d","sig","d_smaller_concord")
names(bl)=names_bl
bl=bl[bl$good_trip==TRUE,]
hist(bl$trl,nclass=1000)        

bl$Order="Odonata"

bl$Suborder=ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Zygoptera),1,all) ,"Zygoptera",
                      ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Epiprocta),1,all) ,"Epiprocta","RANDOM"))
                     

bl$Superfamily_between=ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Anisozygoptera),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Anisoptera),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Zygoptera),1,any) ,"1",
                           ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Aeshnidae),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",c34),1,any),"2",
                           ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Gomphidae_Petaluridae),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",c4),1,any),"3",       
                           ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Cordulegastroidea),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Libelluloidea),1,any),"4",      
                           ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Lestoidea),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",c86),1,any),"5",
                           ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Platystictidae),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",c6),1,any),"8",
                           ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Coenagrionoidea),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Calopterygoidea),1,any),"6",
                           ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Anisozygoptera),1,any) & apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Anisoptera),1,any),"7","RANDOM"))))))))      
        

bl$Superfamily_within=ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Aeshnidae),1,all),"Aeshnoidea(Aeshnidae)",
                       ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Gomphidae_Petaluridae),1,all),"Aeshnoidea(Gomphidae+Petaluridae)",       
                       ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(bl[,c("P1","P2","P3")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM")))))))        
#Distance        
       
        

ch2=ggplot(bl[bl$Order!="RANDOM" & bl$trl<15 & bl$sig==T & bl$d_smaller_concord==T,], aes(x=d, fill=topo))+geom_density(alpha=0.5)+scale_fill_manual(values=c("dodgerblue4", "gold","orange"))+facet_wrap(~Order)+xlab("distance")+ theme(legend.position = "none")+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1))+xlim(0,0.2) 
ch3=ggplot(bl[bl$Suborder!="RANDOM" & bl$trl<15 & bl$sig==T & bl$d_smaller_concord==T,], aes(x=(brl1+brl2)/(2*trl), fill=topo))+geom_density(alpha=0.5)+scale_fill_manual(values=c("dodgerblue4", "gold","orange"))+facet_wrap(~Suborder)+xlab("")+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1),legend.position=c(0.5,-0.6),legend.direction="horizontal",legend.background = element_blank())+xlim(0,0.1)
ch4=ggplot(bl[bl$Superfamily_between!="RANDOM" & bl$trl<10 & bl$trl<15 & bl$sig==T & bl$d_smaller_concord==T,], aes(x=d, fill=topo))+geom_density(alpha=0.5)+scale_fill_manual(values=c("dodgerblue4", "gold","orange"))+facet_wrap(~Superfamily_between)+xlab("")+ theme(legend.position = "none")+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1))+xlim(0,0.15)
ch5=ggplot(bl[bl$Superfamily_within!="RANDOM" & bl$trl<15 & bl$sig==T & bl$d_smaller_concord==T,], aes(x=(brl1+brl2)/(2*trl), fill=topo))+geom_density(alpha=0.5)+scale_fill_manual(values=c("dodgerblue4", "gold","orange"))+facet_wrap(~Superfamily_within)+xlab("")+ theme(legend.position = "none")+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1))+xlim(0,0.05)        
        
        
blaz=bl[bl$Superfamily_between=="1",]
blaz$Superfamily_between="Anisozygoptera"
ch_az=ggplot(blaz[blaz$sig==T & blaz$d_smaller_concord==T & blaz$trl< 15,], aes(x=d, fill=topo))+geom_density(alpha=0.5)+scale_fill_manual(name = "",values=c("dodgerblue4", "gold","orange"))+xlim(0,0.2)+facet_wrap(~Superfamily_between)+xlab("")+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1),legend.position=c(0.8,0.9),legend.direction="vertical",legend.background = element_blank(),legend.key.size = unit(0.5, "cm"),legend.key.width = unit(0.5,"cm"))+ggtitle("C")+xlab("distance")
        
        +geom_vline(data=mu, aes(xintercept=grp.mean, color=topo),linetype="dashed",size=1)+scale_color_manual(name = "",values=c("dodgerblue4", "gold","orange"))   

        
        
lay = rbind(c(1,1,1,1,1),
            c(2,2,NA,NA,NA))
            

f=arrangeGrob(ch1,ch_az,layout_matrix = lay)
ggsave("Fig5_source.png",f, width = 7, height = 6)        
        
lay = rbind(c(1,1,2,2),
            c(3,3,3,3),
           c(3,3,3,3),
           c(3,3,3,3),
           c(4,4,4,4))        

f=arrangeGrob(ch2,ch3,ch4,ch5,layout_matrix = lay)
ggsave("Fig5Suppl.png",f, width = 7, height = 7)          
        
        
        
ggplot(bl[1:10000,],aes((brl1/2)+(brl2/2),brl_int)) + geom_point(size=0.4,stroke=0)+geom_smooth(method = "auto", size = 0.5,color="red")
        
        
        
################################################################ QuIBL/Chi-square ########################################################################

tt=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")

Aeshnoidea=c(extract.clade(tt,137)$tip.label,extract.clade(tt,144)$tip.label)
Calopterygoidea=extract.clade(tt,91)$tip.label
Coenagrionoidea=extract.clade(tt,109)$tip.label
Lestoidea=extract.clade(tt,131)$tip.label
Cordulegastroidea=extract.clade(tt,151)$tip.label
Libelluloidea=extract.clade(tt,156)$tip.label

Zygoptera=extract.clade(tt,88)$tip.label
Anisozygoptera="Epiophlebia_superstes"
Anisoptera=extract.clade(tt,136)$tip.label

Aeshnidae=extract.clade(tt,137)$tip.label
Gomphidae_Petaluridae=extract.clade(tt,144)$tip.label
Libellulidae=extract.clade(tt,161)$tip.label



names_v=c("suborder","id","triplet","outgroup","C1","C2","mixprop1","mixprop2","lambda2Dist","lambda1Dist","BIC2Dist","BIC1Dist","count")
total=read.csv("quibl_results.txt",stringsAsFactors=FALSE,header=FALSE)
names(total)=names_v
P1=gsub(" ","_",paste(unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",1)),unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",2))))
P2=gsub(" ","_",paste(unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",3)),unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",4))))
P3=gsub(" ","_",paste(unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",5)),unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",6))))
total$P1=P1
total$P2=P2
total$P3=P3



#Chisq test for extreme ILS
pvals=c()
common=c()
i=0
for (trip in unique(as.character(total$triplet)))
{
    i=i+1
    progress(i,length(unique(as.character(total$triplet))))
    gcount=total[total$triplet==trip,"count"]+1
    p_v_all=chisq.test(gcount)$p.value
    pos=rank(gcount,ties.method ="random")
    cl=c("discord2","discord1","common")
    common=c(common,cl[pos])
    p_v_d=chisq.test(gcount[pos!=3])$p.value
    pvals=c(pvals,c(NA,p_v_d,p_v_all)[pos])
}

#FDR correction 
total$common=common
total$Q=pvals
total[total$common=="common","Q"]=p.adjust(total[total$common=="common","Q"],method="fdr")
total[total$common=="discord1","Q"]=p.adjust(total[total$common=="discord1","Q"],method="fdr")



#delta BIC
total$BICdiff = total$BIC2-total$BIC1
total[is.na(total$type),"BICdiff"]=0
#Taxa assignment
total$Order="Odonata"

total$P1subo=ifelse(total$P1 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P1 %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$P3subo=ifelse(total$P3 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P3 %in% Anisoptera,"Anisoptera","Anisozygoptera"))

total$P2subo=ifelse(total$P2 %in% Zygoptera,"Zygoptera",
                    ifelse(total$P2 %in% Anisoptera,"Anisoptera","Anisozygoptera"))



total$Suborder=ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Zygoptera),1,all),"Zygoptera",
                      ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Anisoptera),1,all),"Anisoptera",
                      ifelse(total$suborder=="Anisozygoptera" & !apply(t(apply(total[,c("P1subo","P2subo","P3subo")],1,duplicated)),1,any) & total$outgroup!="Epiophlebia_superstes","Anisozygoptera","RANDOM")))

total$Superfamily=ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM")))))

total$focalclade=ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Aeshnidae),1,all),"Aeshnidae",
                      ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Gomphidae_Petaluridae),1,all),"Gomphidae+Petaluridae",
                      ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Libellulidae),1,all),"Libellulidae","RANDOM")))

total$Order=ifelse(total$Suborder!="RANDOM","Odonata","RANDOM")

total$type=ifelse(total$common=="common" & total$BICdiff< -30 ,"Concordant",ifelse(total$BICdiff< -30 & total$common!="common","Introgression+ILS",ifelse(total$common=="common" & total$BICdiff> -30,"Extreme ILS","ILS")))


total$Qtype=ifelse(total$common=="common" & total$Q>10^-6,"Extreme ILS",ifelse(total$common=="discord1" & total$Q>10^-6,"ILS",ifelse(total$common=="discord1" & total$Q<10^-6,"Introgression+ILS","none")))
for (trip in unique(as.character(total$triplet)))
{
    if(total[total$triplet==trip & total$common=="common","Qtype"]=="Extreme ILS")
    {
        total[total$triplet==trip & total$common=="discord1","Qtype"]="none"     
    }    
}



#Mixprop distributions 
total_mixprop=melt(total[total$type=="Introgression+ILS" ,c("mixprop2","Order","Suborder","Superfamily","focalclade")],value.name="taxon",id=c("mixprop2"))
total_mixprop=total_mixprop[total_mixprop$taxon!="RANDOM",]
total_mixprop$variable=replace(as.character(total_mixprop$variable),as.character(total_mixprop$variable)=="focalclade","Focal clade")
total_mixprop$variable_f=factor(total_mixprop$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))
total_mixprop=total_mixprop[complete.cases(total_mixprop), ]

###############Plotting
m1=ggplot(total_mixprop, aes(x=taxon, y=log(mixprop2)))+geom_violin(fill="rosybrown2")+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1))+ylab(expression(log(pi[2])))+xlab("")+ggtitle("A")

#QuIBL proportions 
total_p=melt(total[ ,c("type","Order","Suborder","Superfamily","focalclade")],id="type",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="focalclade","Focal clade")
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))
total_p=total_p[complete.cases(total_p), ]

###############Plotting
m2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=type))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5),size=3)+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1),legend.position=c(0.5,1.4),legend.direction="horizontal",legend.background = element_blank())+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("wheat","gold","grey50","rosybrown2"),name="")+ggtitle("B")


#Chi-square proportions
total_p=melt(total[ ,c("Qtype","Order","Suborder","Superfamily","focalclade")],id="Qtype",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM" & total_p$Qtype!="none",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="focalclade","Focal clade")
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))
total_p=total_p[complete.cases(total_p), ]

mq2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=Qtype))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5),size=3)+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1),legend.position=c(0.5,1.4),legend.direction="horizontal",legend.background = element_blank())+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("gold","grey50","wheat"),name="")+ggtitle("C")




quartz(width=7,height=2.4) 
grid.arrange(mq2,nrow=1)
quartz.save("Fig5.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("Fig5.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)

quartz(width=7,height=4.8) 
grid.arrange(m1,m2,nrow=2)
quartz.save("quibl_mixprop.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("quibl_mixprop.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)

#Introgressing pair
get_intropair_quibl=function(m)
{
    pair_v=c()
    for (i in 1:nrow(m))        
    {
        pair_intro=sort(as.character(total[i,c("P1","P2","P3")][!total[i,c("P1","P2","P3")] %in% total[i,"outgroup"]]))
        pair_v=rbind(pair_v,pair_intro)
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1","i2")
    return(pair_v)
}    

intropairs=get_intropair_quibl(total)
i1i2=as.vector(apply(intropairs,1,paste,collapse="_"))
total_quibl=cbind(total,intropairs)    
total_quibl$i1i2=i1i2



################################################################ Additional Figures ########################################################################
############################################################################################################################################################
############################################################################################################################################################

hyde_bar, g2 , mq2

lay = rbind(c(1,1,1),
             c(2,2,2),
           c(3,3,3))

quartz(width=7,height=7) 
grid.arrange(hyde_bar,g2,mq2,layout_matrix = lay)

quartz.save("Fig3.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("Fig3.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)



#Supplement
lay = rbind(c(1,1,1),
             c(2,2,2),
           c(3,4,5),
           c(6,6,6))

quartz(width=7,height=9.5) 
grid.arrange(a1 ,a2 ,a4 ,a5 ,a6 ,g1,layout_matrix = lay)

quartz.save("Suppl3.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("Suppl3.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)





################################################################ Analysis ##################################################################################
############################################################################################################################################################
############################################################################################################################################################

################################################################ Venn Overlap ##############################################################################

#Odonata
hyde_u=unique(total_hyde[total_hyde$introgressionid=="Introgression","i1i2"])
dfoil_u=unique(total_dfoil[total_dfoil$introgressionid=="Inter-group" & total_dfoil$i1i2!="none_none" ,"i1i2"])
quibl_u=unique(total_quibl[total_quibl$Qtype=="Introgression+ILS" ,"i1i2"])
venn.diagram(x = list(quibl_u,dfoil_u,hyde_u),category.names = c("Chi-square" , "Dfoil","HyDe"),filename = "Odonatavenn.png",output=T,resolution =300)

#Suborder
for (i in unique(total_hyde$Suborder))
{
    hyde_u=unique(total_hyde[total_hyde$introgressionid=="Introgression" & total_hyde$Suborder==i,"i1i2"])
    dfoil_u=unique(total_dfoil[total_dfoil$introgressionid=="Inter-group" & total_dfoil$i1i2!="none_none" & total_dfoil$Suborder==i,"i1i2"])
    quibl_u=unique(total_quibl[total_quibl$Qtype=="Introgression+ILS" & total_quibl$Suborder==i ,"i1i2"])
    venn.diagram(x = list(quibl_u,dfoil_u,hyde_u),category.names = c("Chi-square" , "Dfoil","HyDe"),filename = paste(i,"venn.png",sep=""),output=T,resolution =300)
}    

#Superfamilies
for (i in unique(total_hyde$Superfamily))
{
    hyde_u=unique(total_hyde[total_hyde$introgressionid=="Introgression" & total_hyde$Superfamily==i,"i1i2"])
    dfoil_u=unique(total_dfoil[total_dfoil$introgressionid=="Inter-group" & total_dfoil$i1i2!="none_none" & total_dfoil$Superfamily==i,"i1i2"])
    quibl_u=unique(total_quibl[total_quibl$Qtype=="Introgression+ILS" & total_quibl$Superfamily==i ,"i1i2"])
    venn.diagram(x = list(quibl_u,dfoil_u,hyde_u),category.names = c("Chi-square" , "Dfoil","HyDe"),filename = paste(i,"venn.png",sep=""),output=T,resolution =300)
}    



#Focal clades
for (i in unique(total_hyde$focalclade))
{
    hyde_u=unique(total_hyde[total_hyde$introgressionid=="Introgression" & total_hyde$focalclade==i,"i1i2"])
    dfoil_u=unique(total_dfoil[total_dfoil$introgressionid=="Inter-group" & total_dfoil$i1i2!="none_none" & total_dfoil$focalclade==i,"i1i2"])
    quibl_u=unique(total_quibl[total_quibl$Qtype=="Introgression+ILS" & total_quibl$focalclade==i ,"i1i2"])
    venn.diagram(x = list(quibl_u,dfoil_u,hyde_u),category.names = c("Chi-square" , "Dfoil","HyDe"),filename = paste(i,"venn.png",sep=""),output=T,resolution =300)
}    

################################################################ TMRCA of introgressing species pairs ##########################################################

###Find Tmrca
findTMRCA=function(table_introg_pairs,dated_tree_path,mcmc_table)
{
    phy_mcmc=readMCMCtree(dated_tree_path)
    phy=phy_mcmc$apePhy
    tree_h=nodeheight(phy, node=85)
    age_v=c()
    node_v=c()
    age_p_v=c()
    total_b=table_introg_pairs
    for (i in 1:nrow(total_b))
    {
        progress(i,nrow(total_b))
        age_pair=tree_h-findMRCA(phy,c(total_b[i,"i1"],total_b[i,"i2"]),type="height")
        node_pair=findMRCA(phy,c(total_b[i,"i1"],total_b[i,"i2"]),type="node")
        age_p_pair=sample(mcmc_table[,as.character(node_pair)],10)
        age_v=c(age_v,age_pair)
        node_v=c(node_v,node_pair)
        age_p_v=rbind(age_p_v,age_p_pair)
    }
    age_p_v=data.frame(age_p_v)
    names(age_p_v)=paste("pp",1:10,sep="")
    return(cbind(data.frame(total_b,TMRCA=age_v,nodeN=node_v),age_p_v))
}    

mcmc_age=read.table("mcmc4long.txt",header=T)
mcmc_age$Gen=NULL
mcmc_age$mu=NULL
mcmc_age$sigma2=NULL
mcmc_age$lnL=NULL
names(mcmc_age)=86:(86-1+ncol(mcmc_age))

hyde_nd=total_hyde[!duplicated(total_hyde$i1i2),]
dfoil_nd=total_dfoil[!duplicated(total_dfoil$i1i2),]
quibl_nd=total_quibl[!duplicated(total_quibl$i1i2),]



t_hyde=findTMRCA(total_hyde,"FigTree4long.tre",mcmc_age)
t_dfoil=findTMRCA(total_dfoil,"FigTree4long.tre",mcmc_age)
t_quibl=findTMRCA(total_quibl,"FigTree4long.tre",mcmc_age)




#Time-introgression plot for the entire tree 
all_a=unlist(t_hyde[,paste("pp",1:10,sep="")])
sig_a=unlist(t_hyde[t_hyde$introgressionid=="Introgression" ,paste("pp",1:10,sep="")])
d_a=data.frame(Age=c(all_a,sig_a),Distribution=c(rep("All",length(all_a)),rep("Sig.",length(sig_a))))
quartz(width=6.5, height=4.3)
ggplot(d_a, aes(x=Age, fill=Distribution))+geom_density(alpha=1,position = "stack")+scale_fill_manual(values=c("dodgerblue4", "gold"))+ geom_vline(xintercept = unique(sig_a),linetype="dashed", color = "red", size=1)




################################################################ Statisitcal tests #############################################################################


###Dfoil tests
get_wilx_dfoil=function(m,stats)
{
    m_out=c()
    for(i in c("Suborder","Superfamily","focalclade"))
    {
       
        for (cl in unique(m[,i]))
        {
            
            my_mean=mean(unlist(m[m[,i]==cl,stats]))
            my_pt=wilcox.test(unlist(m[m[,i]==cl,stats]),unlist(m[m[,i]!=cl,stats]),alternative="two.sided")$p.value
            my_t=t.test(m[m[,i]==cl,stats],mu=0)$p.value
            m_out=rbind(m_out,c(i,cl,my_mean,my_pt,my_t))
        }
            
    }
    m_out=as.data.frame(m_out)
    names(m_out)=c("rank","taxon","average","p_two","p_t")
    m_out=m_out[m_out$taxon!="RANDOM",]
    return(m_out)
        
}

get_wilx_dfoil(total_dfoil,c("DFO_stat","DIL_stat","DFI_stat","DOL_stat"))


#Excess of introgression cases for Anisozygoptera
fisher.test(matrix(c(16096,91648,5896,32456),ncol=2))



################################################################ Branch Length Test ##########################################################
#Rscript blt_anisozygoptera.r BUSCO50_dna_pasta_nopart_iqtree_root.tre BUSCO50_dna_pasta_iqtree_all "Anizo"
#Rscript blt_random.r BUSCO50_dna_pasta_nopart_iqtree_root.tre BUSCO50_dna_pasta_iqtree_all "Random"
#Rscript blt_similar_age.r BUSCO50_dna_pasta_nopart_iqtree_root.tre BUSCO50_dna_pasta_iqtree_all ZYG 131 89 
#Rscript blt_similar_age.r BUSCO50_dna_pasta_nopart_iqtree_root.tre BUSCO50_dna_pasta_iqtree_all ANISO 137 143 


tt=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")
Zygoptera=extract.clade(tt,88)$tip.label
Anisoptera=extract.clade(tt,136)$tip.label
bld=read.table("brls_odo.txt",header=T)
bld$clade=rep(c("Anisozygoptera","Random","Anisoptera","Zygoptera"),c(1185005,1896785,1853621,2900520))

ggplot(bld, aes(x=((brl1/trl)+(brl2/trl))/2, fill=topo))+geom_density(alpha=0.5)+scale_fill_manual(values=c("dodgerblue4", "gold","orange"))+xlim(0,0.15)+facet_wrap(~clade)+xlab("distance s1 s2")
ggplot(bld, aes(x=brl_int/trl, fill=topo))+geom_density(alpha=0.5)+scale_fill_manual(values=c("dodgerblue4", "gold","orange"))+xlim(0,0.15)+facet_wrap(~clade)+xlab("internal branch length")

#bld$topo=ifelse(bld$root_tip %in% Zygoptera, "concordant",ifelse(bld$root_tip %in% Anisoptera,"discord1","discord2"))



ggplot(bld, aes(x=((brl1/trl)+(brl2/trl))/2, fill=topo))+geom_density(alpha=1,position = "stack")+scale_fill_manual(values=c("dodgerblue4", "gold","orange"))+xlim(0,0.1)


names_vb=c("clade","P1out","P2out","P3out","CountP1","CountP2","CountP3","PvalueChi","meanT_concord","meanT_discord1","meanT_discord2","PvalueWCOMC1","PvalueWCOMC2","PvalueWC1C2")
total_b=read.csv("Anizo_blt.txt",stringsAsFactors=FALSE,header=F)
names(total_b)=names_vb
total_b$PvalueChi=p.adjust(total_b$PvalueChi,method="fdr")
total_b$PvalueWCOMC1=p.adjust(total_b$PvalueWCOMC1,method="fdr") 
total_b$PvalueWCOMC2=p.adjust(total_b$PvalueWCOMC2,method="fdr") 
total_b$PvalueWC1C2=p.adjust(total_b$PvalueWC1C2,method="fdr") 

b_wilcox=get_intropair_blt_wilcox(total_b)
b_chisq=get_intropair_blt_chisq(total_b)

m_ch=b_chisq
m_wilx=b_wilcox
m_overlap=cbind(b_chisq[,1:2],b_chisq$pass=="TRUE" & b_wilcox$pass=="TRUE")
names(b_wilcox)=c("i1_wilx","i2_wilx","pass_wilx","totalbl_pass")
names(b_chisq)=c("i1_chi","i2_chi","pass_chi")
names(m_overlap)=c("i1","i2","pass")




