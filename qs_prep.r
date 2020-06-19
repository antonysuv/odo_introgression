library("ggplot2")
library("ape")
#Quartet Sampling analysis 
files = list.files(".", pattern="*.scores", full.names=TRUE, recursive=FALSE)
all_scores=rbind()
for (f in files)
{
    tt=read.csv(f)
    if (tt[1,1]=="QS1")
    {
        tt=tt[3:83,]
    }
    else
    {
        tt=tt[1:81,]
    }
    tt=t(tt)[2:5,]
    tt=apply(tt,2,as.numeric)
    row.names(tt)=c("freq0","qc","qd","qi")
    all_scores=rbind(all_scores,tt)
}    

#Aniso
nodes=as.character(c(3,4,5,6,25,37,38,39,62,63,71))
n_names=c("Anisoptera+\nAnisozygoptera","Aeshnidae+\nExophytica","Cavilabiata+\n(Gomphidae+\nPetaluridae)","Cordulegastroidea+\nLibelluloidea","Gomphidae+\nPetaluridae",
         "Lestoidea+\n(Platystictoidea+\nCalopterigoidea+\nCoenagrionoidea)","Platystictoidea+\n(Calopterigoidea+\nCoenagrionoidea)","Calopterigoidea+\nCoenagrionoidea","(A+B)+\n(C+D)","A+B","C+D")
sc=matrix(all_scores[,nodes])
stat=rep(row.names(all_scores[,nodes]))
n=as.vector(rep(n_names,each=nrow(all_scores[,nodes])))
zz=data.frame(score=sc,statistic=stat,node=n)


p = ggplot(zz, aes(x=node, y=score,fill=statistic)) + geom_boxplot() + 
     scale_x_discrete(limits=n_names) +
     scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","white"),name = "Score", labels = c("Frequency", "QC", "QD","QI"))+
     theme_classic()+
     labs(x="Node", y = "Score value")+
     theme(axis.text.x = element_text(angle = 360,face = "bold")) 
     xlab(1:5)



phy=read.tree("/Users/Anton/Downloads/RUN_odoqs.labeled.tre")
phy=drop.tip(phy,"Ephemera_danica")
phy$node.label=1:83
brtimes=sort(branching.times(phy),decreasing=T)
n_order=names(brtimes)[3:83]
qd=all_scores[row.names(all_scores)=="qd",n_order]
f=all_scores[row.names(all_scores)=="freq0",n_order]
qd_min=apply(qd,2,min,na.rm = T)
qd_max=apply(qd,2,max,na.rm = T)

qd_traj=rbind()
for (i in 1:1000)
{
    mm=apply(rbind(apply(qd,2,function(x){sample(x[!is.na(x)], size = 1)}),apply(qd,2,function(x){sample(x[!is.na(x)], size = 1)})),2,sort)
    y=runif(81,mm[1,],mm[2,])
    qd_traj=rbind(qd_traj,y)
}   

plot(brtimes[3:83],qd_traj[1,],type="l")
apply(qd_traj,1,lines,type="l",x=brtimes[3:83],col=adjustcolor( "navy", alpha.f = 0.5))
