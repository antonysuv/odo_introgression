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
