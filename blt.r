library("doSNOW")
library("ape")
cl = makeCluster(30,type = "SOCK") 


library("phangorn")
library("foreach")
#library("doParallel")



args = commandArgs(trailingOnly=TRUE)




test_triplet=function(taxa,gene_trees,clade_name,sp_tree)
{
    
    trl_all=c()
    brls_all1=c()
    brls_all2=c()
    internal_all=c()
    out_all=c()
    root_tip_all=c()
    concord_all=c()
    tripl_all=c()
    ind=0
    sp_trip=keep.tip(sp_tree,taxa)
    sp_out=drop.tip(sp_trip,extract.clade(sp_trip,node=5)$tip.label)$tip.label
    for (tre in gene_trees)
    {
        
        ind=ind+1
        if(any(c("Ephemera_danica","Isonychia_kiangsinensis") %in% tre$tip.label) & all(taxa %in% tre$tip.label))
        {
            
            tre=root(tre,sample(tre$tip.label[tre$tip.label %in% c("Ephemera_danica","Isonychia_kiangsinensis")],1))
            trl=sum(tre$edge.length)
            tre_trip=keep.tip(tre,taxa)
            concord_t=as.numeric(drop.tip(tre_trip,extract.clade(tre_trip,node=5)$tip.label)$tip.label!=sp_out)
            brls=extract.clade(tre_trip,max(tre_trip$edge))$edge.length
            root_tip=tre_trip$tip.label[!tre_trip$tip.label %in% extract.clade(tre_trip,max(tre_trip$edge))$tip.label]
            trl_all=c(trl_all,trl)
            outl=tre_trip$edge.length[which(tre_trip$edge[,1]==4 & (tre_trip$edge[,2]==1 | tre_trip$edge[,2]==2 | tre_trip$edge[,2]==3))]
            out_all=c(out_all,outl)
            internall=tre_trip$edge.length[which(tre_trip$edge[,1]==4 & tre_trip$edge[,2]==5)]
            internal_all=c(internal_all,internall)
            brls_all1=c(brls_all1,brls[1])
            brls_all2=c(brls_all2,brls[2])
            tripls=mean(c(mean(c(brls[1],brls[2]))+internall,outl))
            tripl_all=c(tripl_all,tripls)
            root_tip_all=c(root_tip_all,root_tip)
            concord_all=c(concord_all,concord_t)
            
        }    
    } 
    counts=table(root_tip_all)
    con=names(which.max(counts))
    dis=names(which.min(counts))
    #0=concord 1=discord1 2=discord2 
    m=data.frame(clade_name,P1=tre_trip$tip.label[1],P2=tre_trip$tip.label[2],P3=tre_trip$tip.label[3],brl1=brls_all1,brl2=brls_all2,trl=trl_all,brl_out=out_all,brl_int=internal_all,root_tip=root_tip_all,topo=ifelse(root_tip_all %in% con,"concord",ifelse(root_tip_all %in% dis,"discord1","discord2")),topo_con_sp=concord_all,good_trip=con==sp_out,tripl=tripl_all)
    if(!any(table(m$root_tip)==0) & length(table(m$root_tip))==3)
    {
        counts=table(m$root_tip)
        com=names(which.max(counts))
        congruent=sp_out==com
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
        v_out=c(clade_name,names(counts),counts,congruent,chi,mean(ccom),mean(c1),mean(c2),w_testc1,w_testc2,w_test)
        m$sig=chi<0.05 & w_test<0.05
        write.table(m,paste(c(taxa,"csv"),collapse="."),quote=F,row.names=F,col.names=F)
        return(as.vector(v_out))
    }

}    


getstats_triplets=function(gene_trees,taxa_combn,clade_name,sp_tree)
{
       
    tt=sp_tree
    print(ncol(taxa_combn))
    pb=txtProgressBar(0,ncol(taxa_combn),style=3)
    progress=function(n){
    setTxtProgressBar(pb,n)
    }
    opts=list(progress=progress)
    
    
    out_t=foreach(i=1:ncol(taxa_combn),.combine='rbind',.options.snow=opts) %dopar% 
    {
       
        triplet=taxa_combn[,i]
        stats=test_triplet(triplet,gene_trees,clade_name,tt)
        return(stats)
        
    }
                        names_vb=c("clade","P1out","P2out","P3out","CountP1","CountP2","CountP3","congruent","PvalueChi","meanT_concord","meanT_discord1","meanT_discord2","PvalueWCOMC1","PvalueWCOMC2","PvalueWC1C2")
    out_t=data.frame(out_t)
    names(out_t)=names_vb
    write.table(as.data.frame(out_t),paste(clade_name,".out",sep=""),quote = F, row.names = F, col.names = F,sep=",")
    
}    

clusterExport(cl, c("root","keep.tip","extract.clade","setTxtProgressBar","test_triplet","drop.tip"))
registerDoSNOW(cl)


tt=read.tree("../BUSCO50_dna_pasta_nopart_iqtree_root.tre")
g_trees=read.tree("../BUSCO50_dna_pasta_iqtree_all")



####Anisozygoptera
clade1=extract.clade(tt,as.numeric(88))$tip.label
clade2=extract.clade(tt,as.numeric(136))$tip.label
clade3="Epiophlebia_superstes"
taxa_combn_az=t(expand.grid(clade1,clade2,clade3))
getstats_triplets(g_trees,taxa_combn_az,"Anisozygoptera",tt)
####Zygoptera
clade_z=tt$tip.label[1:48]
taxa_combn_z=combn(clade_z,3)
getstats_triplets(g_trees,taxa_combn_z,"Zygoptera",tt)
####Epiprocta
clade_a=tt$tip.label[49:83]
taxa_combn_a=combn(clade_a,3)
getstats_triplets(g_trees,taxa_combn_a,"Epiprocta",tt)


stopCluster(cl)