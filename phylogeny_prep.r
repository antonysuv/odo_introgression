library(MCMCtreeR)

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

MCMC.tree.plot(phy_mcmc,analysis.type = "MCMCtree",MCMC.chain =tt,plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), no.margin = T, cex.labels = 0.01,cex.tips = 0.6,ladderize.tree = F,pos.age=-7,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6)

MCMC.tree.plot(phy_mcmc,analysis.type = "MCMCtree",MCMC.chain =tt,plot.type = "cladogram",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), no.margin = T, cex.labels = 0.001,cex.tips = 0,ladderize.tree = F,pos.age=-7,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6, show.tip.label = F)