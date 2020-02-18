library('ape')
library('phytools')
library('phangorn')



#ML Phylonet
phylonet_chunk_topo=function(sp_tree,gene_trees,sp_list,n_retic)
{
    sp_phy=read.tree(sp_tree)
    sp_phy$node.label=NULL
    sp_phy_sub=keep.tip(sp_phy,sp_list[sp_list %in% sp_phy$tip.label])
    phy=read.tree(gene_trees)
    rooted_phy=c()
    out_v=c("Isonychia_kiangsinensis","Ephemera_danica")
    for (gt in phy)
    {
        if (any(sp_list %in% gt$tip.label))
        {
             gt_sub=keep.tip(gt,sp_list[sp_list %in% gt$tip.label])
        } 
        
        if (any(out_v %in% gt_sub$tip.label))
        {
            out_sp=out_v[which(out_v %in% gt_sub$tip.label)[1]]
            gt_r=root(gt_sub,out_sp,resolve.root=T)
            gt_r$node.label[gt_r$node.label=="Root"]=100
            gt_r$node.label[gt_r$node.label==""]=1
            gt_r$node.label=as.numeric(gt_r$node.label)/100
            rooted_phy=c(rooted_phy,write.tree(gt_r))
                
        } 
        
       
    } 
    phy=rooted_phy
    for (i in 1:n_retic)
    {    
        f_n=paste("phylonet_genes_",i,"ret.nex",sep="")
        write("#NEXUS\n\nBEGIN TREES;",f_n)
        write(paste("Tree fixtr = ",write.tree(sp_phy_sub)),f_n,append=T)
        d=data.frame(rep("Tree",length(phy)),paste("gt",1:length(phy),"=",sep=""),phy)
        write.table(d,f_n,quote = F,row.names = F, col.names=F,append=T)
        write("END;\n\nBEGIN PHYLONET;",f_n,append=T)
        write(paste("InferNetwork_ML (all)",i,"-s fixtr -fs -di -pl 15 -x 1;","\nEND;"),f_n,append=T)
    }    
} 


anax=c('Aeshna_palmata','Anax_junius','Anax_walsinghami','Anax_parthenope','Gynacantha_tibiata','Austroaeschna_subapicalis','Telephlebia_godeffroyi',"Isonychia_kiangsinensis","Ephemera_danica") 
gompeta=c('Asiagomphus_melaenops','Phanogomphus_spicatus','Stylurus_spiniceps','Leptogomphus_perforatus','Ictinogomphus_pertinax','Phenes_raptor','Tanypteryx_pryeri','Gomphomacromia_paradoxa',"Isonychia_kiangsinensis","Ephemera_danica")

phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",anax,1)
phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",gompeta,1)