library(MCMCtreeR)
library(ape)
library(laser)
library(geiger)
library(TreeSim)
library(phytools)
library(ggplot2)
library(reshape2)
library(gridExtra)
source("http://ib.berkeley.edu/courses/ib200b/labs/lab12/rbdtree.n3.R")

phy_mcmc=readMCMCtree("/Users/Anton/Downloads/all_dated_trees/run4long.tree.tre")
phy=phy_mcmc$apePhy
phy=rotate(phy,86)
phy=rotate(phy,87)
phy=rotate(phy,135)
phy=rotate(phy,136)
phy=rotate(phy,90)
phy=rotate(phy,91)
phy=rotate(phy,143)
phy=rotate(phy,110)
phy=rotate(phy,101)
phy=rotate(phy,92)
phy=rotate(phy,150)
phy=rotate(phy,156)
phy=rotate(phy,158)
phy=rotate(phy,159)
phy=rotate(phy,161)
phy=rotate(phy,163)
phy=rotate(phy,167)
phy=rotate(phy,118)
phy_mcmc$apePhy=phy

tt=read.table("/Users/Anton/Downloads/mcmc.txt",header=T)
tt=tt[,1:85]

# Fig 1 Main dated tree
MCMC.tree.plot(phy_mcmc,analysis.type = "MCMCtree",MCMC.chain =tt,plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), no.margin = T, cex.labels = 0.01,cex.tips = 0.6,ladderize.tree = F,pos.age=-7,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6)

MCMC.tree.plot(phy_mcmc,analysis.type = "MCMCtree",MCMC.chain =tt,plot.type = "cladogram", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), no.margin = T, cex.labels = 0.001,cex.tips = 0,ladderize.tree = F,pos.age=-300,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6, show.tip.label = F)

# Fig 3 LTT
tt=extract.clade(phy,87)
yule_odo=yule(tt)
btimes=sort(branching.times(tt),decreasing=TRUE)
#83 is the number of tips
tt_null=replicate(10000,rbdtree.n(83,btimes[1],yule_odo$lambda,0),simplify=FALSE)
mycol=rgb(190, 190, 190, max = 255, alpha = 10)
#Read fossil data
fos=read.csv("/Users/Anton/Desktop/git_repos/macrology/data/pbdb2019_datasets/odonatoptera_pbdb2019.csv" ,sep=",")
fos=fos[fos$accepted_rank=="species",]
fos_dens=apply(fos[,c("max_ma","min_ma")],1,mean)

#Plotting
par(mar = c(5,5,2,5))
ltt.plot(tt_null[[1]],col=mycol,log="y",xlim=c(-320,0),ylim=c(1,90),ylab="Extant lineages",xlab="Time (Myr ago)")
for (i in 2:10000)
{
        ltt.lines(tt_null[[i]],col=mycol)
        #gam=c(gam,ltt(tt_null[[i]],log=F,plot=F)$gamma)
}       
ltt.lines(tt,col="black",lwd=2)
lines(-131.9361,10,type="p",pch=16,col="red")

par(new = T)
hist(-(fos_dens),xlim=c(-320,0),axes=F, xlab=NA, ylab=NA,border=F,nclass=30,col=rgb(138, 163, 255, max = 255, alpha = 90),main="",ylim=c(0,250))
axis(side = 4)
mtext(side = 4, line = 3, 'Fossil abundance')
#P-Tr
abline(v=-252,lty=2)
#Tr-J
abline(v=-201.3,lty=2)
#K-Pg
abline(v=-66,lty=2)

mtext(side = 3, "P-Tr",at=-252)
mtext(side = 3, "Tr-J",at=-201.3)
mtext(side = 3, "K-Pg",at=-66)
rect(xleft=-145, ybottom=-10, xright=-66, ytop=0,col=rgb(129,196,85, max = 255))
text(x=-105.5, y = -5,"Cretaceous",cex=0.9)
rect(xleft=-145, ybottom=0, xright=-66, ytop=270,col=rgb(129,196,85, max = 255,alpha=40),lty=0)




