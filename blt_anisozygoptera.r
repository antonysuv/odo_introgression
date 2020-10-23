library("ape")
library("phangorn")
library("svMisc")

args = commandArgs(trailingOnly=TRUE)

test_triplet=function(taxa,gene_trees,clade_name)
{
    
    trl_all=c()
    brls_all1=c()
    brls_all2=c()
    internal_all=c()
    out_all=c()
    root_tip_all=c()
    ind=0
    for (tre in gene_trees)
    {
        
        ind=ind+1
        if(any(c("Ephemera_danica","Isonychia_kiangsinensis") %in% tre$tip.label) & all(taxa %in% tre$tip.label))
        {
            
            tre=root(tre,sample(tre$tip.label[tre$tip.label %in% c("Ephemera_danica","Isonychia_kiangsinensis")],1))
            trl=sum(tre$edge.length)
            tre_trip=keep.tip(tre,taxa)
            brls=extract.clade(tre_trip,max(tre_trip$edge))$edge.length
            root_tip=tre_trip$tip.label[!tre_trip$tip.label %in% extract.clade(tre_trip,max(tre_trip$edge))$tip.label]
            trl_all=c(trl_all,trl)
            out_all=c(out_all,tre_trip$edge.length[1])
            internal_all=c(internal_all,tre_trip$edge.length[2])
            brls_all1=c(brls_all1,brls[1])
            brls_all2=c(brls_all2,brls[2])
            root_tip_all=c(root_tip_all,root_tip)
            
        }    
    }    
    counts=table(root_tip_all)
    con=names(which.max(counts))
    dis=names(which.min(counts))
    m=data.frame(brl1=brls_all1,brl2=brls_all2,trl=trl_all,brl_out=out_all,brl_int=internal_all,root_tip=root_tip_all,topo=ifelse(root_tip_all %in% con,"concord",ifelse(root_tip_all %in% dis,"discord2","discord1")))
    write.table(m,paste(c(taxa,"csv"),collapse="."),quote=F,row.names=F)
    if(!any(table(m$root_tip)==0) & length(table(m$root_tip))==3)
    {
        counts=table(m$root_tip)
        com=names(which.max(counts))
        not_com=names(counts)[!names(counts) %in% com]
        not_com_c=counts[!names(counts) %in% com]
        m$common=ifelse(m$root_tip==com,"TRUE","FALSE")
        m$proxy_t=(m$brl1+m$brl2)/m$trl
        ccom=m[m$common==TRUE,"proxy_t"]
        c1=m[m$root_tip==not_com[1],"proxy_t"]
        c2=m[m$root_tip==not_com[2],"proxy_t"]
        w_testc1=wilcox.test(ccom,c1)$p.value
        w_testc2=wilcox.test(ccom,c2)$p.value
        w_test=wilcox.test(c1,c2)$p.value
        chi=chisq.test(not_com_c)$p.value
        v_out=c(clade_name,names(counts),counts,chi,mean(ccom),mean(c1),mean(c2),w_testc1,w_testc2,w_test)
        return(as.vector(v_out))
    } 
  

}    


getstats_triplets=function(phy,gene_trees,clade_name)
{
   
    clade1=extract.clade(tt,as.numeric(88))$tip.label
    clade2=extract.clade(tt,as.numeric(136))$tip.label
    clade3="Epiophlebia_superstes"
    taxa_combn=t(expand.grid(clade1,clade2,clade3))
    
    out_t=c()
    for (i in 1:ncol(taxa_combn))
    {
        progress(i,ncol(taxa_combn))
        triplet=taxa_combn[,i]
        stats=test_triplet(triplet,gene_trees,clade_name)
        out_t=rbind(out_t,stats) 
    }
    write.table(as.data.frame(out_t),clade_name,quote = F, row.names = F, col.names = F,sep=",")
    
}    



tt=read.tree(args[1])
g_tree=read.tree(args[2])
getstats_triplets(tt,g_tree,args[3])
