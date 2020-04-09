library("ape")
library("svMisc")


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



tt=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")
Outgroup=c("Ephemera_danica")
Anisozygoptera=keep.tip(tt,c(extract.clade(tt,131)$tip.label,extract.clade(tt,137)$tip.label,"Epiophlebia_superstes",Outgroup))
Zygoptera=keep.tip(tt,c(extract.clade(tt,88)$tip.label,Outgroup))
Anisoptera=keep.tip(tt,c(extract.clade(tt,136)$tip.label,Outgroup))

dfoil_select(Anisozygoptera,"Anisozygoptera_dfoil")
dfoil_select(Anisoptera,"Anisoptera_dfoil")
dfoil_select(Zygoptera,"Zygoptera_dfoil")

