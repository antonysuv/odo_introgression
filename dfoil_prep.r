library("ape")
library("svMisc")
library('Rtsne')
library('MASS')
library('gridExtra')
library('ggplot2')
library('pals')
library('reshape2')



dfoil_select=function(phy,id)
{
    dfoil_out=c()
    taxa_combn=combn(phy$tip.label,m=5)
    test_phy=read.tree(text="((('A','A'),('A','A')),A);")
    for (i in 1:ncol(taxa_combn))
    {
        progress(i,ncol(taxa_combn))
        sub_phy=keep.tip(phy,taxa_combn[,i])
        sub_phy_dup=sub_phy
        sub_phy_dup$tip.label=rep("A",length(sub_phy_dup$tip.label))
        if (all.equal.phylo(sub_phy_dup,test_phy,use.edge.length=F,use.tip.label=F))
        {
            node_number_order=rev(order(node.depth.edgelength(sub_phy)[(length(sub_phy$tip.label)+1):length(node.depth.edgelength(sub_phy))])+length(sub_phy$tip.label))[1:2]
            taxa_order=c()
            for (n in node_number_order)
            {
               node_taxa=extract.clade(sub_phy,n)$tip.label
               taxa_order=c(taxa_order,node_taxa) 
            }
            dfoil_order=c(taxa_order,sub_phy$tip.label[!sub_phy$tip.label %in% taxa_order])
            dfoil_out=c(dfoil_out,paste(dfoil_order,collapse=","))
            
        }    
        
    }
    write(dfoil_out,id)
}    

dfoil_select_az=function(phy,id)
{
    dfoil_out=c()
    taxa_combn=combn(phy$tip.label,m=5)
    taxa_combn=data.frame(t(taxa_combn))
    taxa_combn=taxa_combn[taxa_combn$X5=="Ephemera_danica",]
    taxa_combn=taxa_combn[apply(taxa_combn=="Epiophlebia_superstes",1,any),]
    write(apply(taxa_combn,1,paste,collapse=","),id)
}    

tt=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")
Outgroup=c("Ephemera_danica")
Zygoptera=keep.tip(tt,c(extract.clade(tt,88)$tip.label,Outgroup))
Anisoptera=keep.tip(tt,c(extract.clade(tt,136)$tip.label,Outgroup))



dfoil_select_az=function(phy,id)
{
    dfoil_out=c()
    taxa_combn=combn(drop.tip(phy,c("Epiophlebia_superstes","Ephemera_danica", "Isonychia_kiangsinensis"))$tip.label,m=3)
    taxa_combn=rbind(taxa_combn,rep("Epiophlebia_superstes",ncol(taxa_combn)))
    taxa_combn=rbind(taxa_combn,rep("Ephemera_danica",ncol(taxa_combn)))
    test_phy=read.tree(text="((('A','A'),('A','A')),A);")
    for (i in 1:ncol(taxa_combn))
    {
        progress(i,ncol(taxa_combn))
        sub_phy=keep.tip(phy,taxa_combn[,i])
        sub_phy_dup=sub_phy
        sub_phy_dup$tip.label=rep("A",length(sub_phy_dup$tip.label))
        if (all.equal.phylo(sub_phy_dup,test_phy,use.edge.length=F,use.tip.label=F))
        {
            node_number_order=rev(order(node.depth.edgelength(sub_phy)[(length(sub_phy$tip.label)+1):length(node.depth.edgelength(sub_phy))])+length(sub_phy$tip.label))[1:2]
            taxa_order=c()
            for (n in node_number_order)
            {
               node_taxa=extract.clade(sub_phy,n)$tip.label
               taxa_order=c(taxa_order,node_taxa) 
            }
            dfoil_order=c(taxa_order,sub_phy$tip.label[!sub_phy$tip.label %in% taxa_order])
            dfoil_out=c(dfoil_out,paste(dfoil_order,collapse=","))
            
        }    
        
    }
    write(dfoil_out,id)
}    



dfoil_select(Anisoptera,"Anisoptera_dfoil")
dfoil_select(Zygoptera,"Zygoptera_dfoil")
dfoil_select_az(tt,"Anisozygoptera_dfoil")

################################################################ RUN DFOIL ######################################################################## 

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
Libellulidae=extract.clade(tt,161)$tip.label



total=read.csv("/Users/Anton/Downloads/dfoil_results.txt",stringsAsFactors=FALSE)
names(total)=names_v
total=total[total$T12<total$T34,]



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
                   ifelse(total$P1 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P1 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P1 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P1 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P1 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))

total$P3sub=ifelse(total$P3 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P3 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P3 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P3 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P3 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P3 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P3 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P3 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P3 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))


total$P2sub=ifelse(total$P2 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P2 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P2 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P2 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P2 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P2 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P2 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P2 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P2 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))


total$P4sub=ifelse(total$P4 %in% Lestoidea,"Lestoidea",
                   ifelse(total$P4 %in% Coenagrionoidea,"Coenagrionoidea",
                   ifelse(total$P4 %in% Calopterygoidea,"Calopterygoidea",
                   ifelse(total$P4 %in% Platystictidae,"Platystictidae",
                   ifelse(total$P4 %in% Aeshnidae,"Aeshnidae",
                   ifelse(total$P4 %in% Gomphidae_Petaluridae,"Gomphidae+Petaluridae",
                   ifelse(total$P4 %in% Cordulegastroidea,"Cordulegastroidea",
                   ifelse(total$P4 %in% Libelluloidea,"Libelluloidea",
                   ifelse(total$P4 %in% Anisozygoptera,"Epiophlebiidae","RANDOM")))))))))


#nodups=total$P1sub!=total$P3sub & total$P1sub!=total$P4sub & total$P2sub!=total$P3sub & total$P2sub!=total$P4sub


total$Order="Odonata"

total$Suborder=ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Zygoptera),1,all) ,"Zygoptera",
                      ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Anisoptera),1,all) ,"Anisoptera",
                      ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Anisozygoptera),1,any) ,"Anisozygoptera","RANDOM")))


total$Superfamily=ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Lestoidea),1,all),"Lestoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Calopterygoidea),1,all),"Calopterygoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Coenagrionoidea),1,all),"Coenagrionoidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Cordulegastroidea),1,all),"Cordulegastroidea",
                       ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Libelluloidea),1,all),"Libelluloidea","RANDOM")))))

total$focalclade=ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Aeshnidae),1,all),"Aeshnidae",
                      ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Gomphidae_Petaluridae),1,all),"Gomphidae+Petaluridae",
                      ifelse(apply(apply(total[,c("P1","P2","P3","P4")],2,"%in%",Libellulidae),1,all),"Libellulidae","RANDOM")))




total$introgressionid=ifelse(total$introgression=="none","None",
                             ifelse(total$introgression=="123" | total$introgression=="124","Ancestral","Inter-group"))



total_ord=melt(total[total$introgressionid!="None",c("DFO_stat","DIL_stat","DFI_stat","DOL_stat","Order","Superfamily","Suborder","focalclade")],measure.vars = c("DFO_stat","DIL_stat","DFI_stat","DOL_stat"),value.name="stat")
total_ord=melt(total_ord,id=c("variable","stat"),value.name="taxon",variable.name="class")
total_ord$class=replace(as.character(total_ord$class),as.character(total_ord$class)=="focalclade","Focal clade")
total_ord$variable_f=factor(total_ord$class, levels=c("Order","Suborder","Superfamily","Focal clade"))
total_ord$variable=ifelse(total_ord$variable=="DFO_stat","D[FO]",ifelse(total_ord$variable=="DIL_stat","D[IL]",ifelse(total_ord$variable=="DFI_stat","D[FI]","D[OL]")))

odo_null=data.frame(variable=rep(c("D[FO]","D[FI]","D[OL]","D[IL]"),3),stat=NA,class=rep(c("Superfamily","Focal clade"),times=c(4,8)),taxon=rep(c("Lestoidea","Aeshnidae","Gomphidae+Petaluridae"),times=c(4,4,4)),variable_f=rep(c("Superfamily","Focal clade"),times=c(4,8)))
total_ord=rbind(total_ord,odo_null)

###############Plotting 
g1=ggplot(total_ord[total_ord$taxon!="RANDOM",], aes(x=taxon, y=stat,fill=variable))+geom_violin(lwd=0.1,bw=0.09)+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="red",position=position_dodge(0.9))+geom_boxplot(width=0.01,outlier.size=-1,position=position_dodge(0.9))+theme(axis.text.x = element_text(size = 8,angle=15,hjust = 1))+ylab("Value")+xlab("")+scale_fill_manual(name="",values=c("black","darkgrey","grey90","white"),labels=expression(D[FI],D[FO],D[IL],D[OL]))+theme(legend.position=c(0.5,1.4),legend.direction="horizontal",legend.background = element_blank())+ggtitle("A") 
                        


total_p=melt(total[,c("introgressionid","Order","Suborder","Superfamily","focalclade")],id="introgressionid",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="focalclade","Focal clade")
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))

g2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=introgressionid))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5),size=3)+theme(axis.text.x = element_text(size = 8, angle=15, hjust = 1),legend.position=c(0.5,1.4),legend.direction="horizontal",legend.background = element_blank())+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("rosybrown2","wheat","grey50"),name="")+ggtitle("B")





#tSNE
total_dr=total[total$introgression!="none",c(8:23)]/apply(total[total$introgression!="none",c(8:23)],1,sum)
tsne_out=Rtsne(total_dr)
sunset=c("#364B9A" ,"#4A7BB7", "#6EA6CD", "#98CAE1" ,"#C2E4EF" ,"#EAECCC", "#FEDA8B" ,"#FDB366", "#F67E4B", "#DD3D2D", "#A50026")
#Main dr figure
dat1=data.frame(tsne_out$Y)
names(dat1)=c("tSNE1","tSNE2")
introg_d=cbind(dat1,introgressionid=total[total$introgression!="none",c("introgressionid")])
#introg_d=introg_d[order(introg_d$introgressionid,decreasing = T),]
dat1=cbind(dat1,total[total$introgression!="none",c("DFO_stat","DIL_stat","DFI_stat","DOL_stat")])
dat1=melt(dat1,id.vars=c("tSNE1","tSNE2"))
dat1$variable=ifelse(dat1$variable=="DFO_stat","D[FO]",ifelse(dat1$variable=="DIL_stat","D[IL]",ifelse(dat1$variable=="DFI_stat","D[FI]","D[OL]")))


t1=ggplot(dat1[dat1$value>-0.5 & dat1$value<0.5,]) + geom_point(aes(tSNE1, tSNE2, color = value),size=0.3,stroke=0)+facet_wrap(~variable,nrow = 1,labeller=label_parsed)+scale_color_gradientn(colors = sunset,name="",guide=guide_colorbar(barwidth = 6,barheight =0.5,label.theme = element_text(size=8)))+theme_classic()+ggtitle("C")+theme(legend.position=c(0.5,1.1) ,legend.direction="horizontal",legend.background = element_blank(),strip.background = element_rect(colour = "white", fill = "white"))



lay = rbind(c(1,1,1,1,1,1),
             c(2,2,2,2,2,2),
           c(3,3,3,3,3,3))


quartz(width=7,height=7.2) 
grid.arrange(g1,g2,t1,nrow=3,layout_matrix = lay)
quartz.save("dfoil_distribution.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("dfoil_distribution.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)








#STATS for Focal clades 

total_Dsig=total[total$introgression!="none",]

get_t=function(m,stats)
{
    m_out=c()
    for(i in c("Order","Suborder","Superfamily","focalclade"))
    {
       
        for (cl in unique(m[,i]))
        {
            
            my_t_two=t.test(m[m[,i]==cl,stats],mu=0)$p.value
            my_t_greater=t.test(m[m[,i]==cl,stats],mu=0,alternative = "greater")$p.value
            my_t_less=t.test(m[m[,i]==cl,stats],mu=0,alternative = "less")$p.value
            
            m_out=rbind(m_out,c(i,cl,my_t_two,my_t_greater,my_t_less))
        }
            
    }
    m_out=as.data.frame(m_out)
    names(m_out)=c("rank","taxon","p_ttest_two","p_ttest_greater","p_ttest_less")
    m_out$psig_two=as.numeric(as.vector(m_out$p_ttest_two))<0.05
    m_out$psig_greater=as.numeric(as.vector(m_out$p_ttest_greater))<0.05
    m_out$psig_less=as.numeric(as.vector(m_out$p_ttest_less))<0.05
    return(m_out[m_out$taxon!="RANDOM",])
        
}






#Total tree test
dat=total[,c("DFO_stat","DIL_stat","DFI_stat","DOL_stat","introgression")]
g1=ggplot(dat)+geom_histogram(aes(DFO_stat),bins=500,col="black")+labs(x = expression(D[FO]),y="N quintets")+ggtitle("A")
g2=ggplot(dat)+geom_histogram(aes(DIL_stat),bins=500,col="black")+labs(x = expression(D[IL]),y="")+ggtitle("B")
g3=ggplot(dat)+geom_histogram(aes(DFI_stat),bins=500,col="black")+labs(x = expression(D[FI]),y="")+ggtitle("C")
g4=ggplot(dat)+geom_histogram(aes(DOL_stat),bins=500,col="black")+labs(x = expression(D[OL]),y="")+ggtitle("D")

dat=dat[total$introgression!="none",c("DFO_stat","DIL_stat","DFI_stat","DOL_stat")]
g5=ggplot(dat)+geom_histogram(aes(DFO_stat),bins=500,col="black")+labs(x = expression(D[FO]),y="N quintets")+ggtitle("E")
g6=ggplot(dat)+geom_histogram(aes(DIL_stat),bins=500,col="black")+labs(x = expression(D[IL]),y="")+ggtitle("F")
g7=ggplot(dat)+geom_histogram(aes(DFI_stat),bins=500,col="black")+labs(x = expression(D[FI]),y="")+ggtitle("G")
g8=ggplot(dat)+geom_histogram(aes(DOL_stat),bins=500,col="black")+labs(x = expression(D[OL]),y="")+ggtitle("H")

grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=4,nrow=2)







dat=dat[dat$Pvalue<10^-6,c("Gamma","D","Pvalue")]
g5=ggplot(dat)+geom_histogram(aes(Gamma),bins=500,col="black")+labs(x = expression(gamma),y="N quartets")+ggtitle("B")
g6=ggplot(dat)+geom_histogram(aes(D),bins=500,col="black")+labs(x = "D",y="N quartets")+ggtitle("")
dat$density = get_density(dat$Gamma, dat$D, n = 100)
g7=ggplot(dat) + geom_point(aes(Gamma, D, color = density),size=0.1)+scale_color_gradientn(colors = viridis(30))+labs(x = expression(gamma),y="D")+ggtitle("")
g8=ggplot(dat) + geom_point(aes(Gamma, D, color = Pvalue),size=0.1)+scale_color_gradientn(colors = viridis(30))+labs(x = expression(gamma),y="D")+ggtitle("")












        
   