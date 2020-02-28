library('ape')
library('phytools')
library('phangorn')



#ML Phylonet
phylonet_chunk_topo=function(sp_tree,gene_trees,sp_list,n_retic)
{
    sp_phy=read.tree(sp_tree)
    sp_phy$node.label=NULL
    out_v=c("Isonychia_kiangsinensis","Ephemera_danica")
    sp_phy_sub=keep.tip(sp_phy,c(sp_list[sp_list %in% sp_phy$tip.label]))
    phy=read.tree(gene_trees)
    rooted_phy=c()
    for (gt in phy)
    {
        if (any(sp_list %in% gt$tip.label) & any(out_v %in% gt$tip.label))
        {
              
            out_sp=out_v[which(out_v %in% gt$tip.label)[1]]
            gt_r=root(gt,out_sp,resolve.root=T)
            gt_sub=keep.tip(gt_r,c(sp_list[sp_list %in% gt$tip.label]))
            gt_sub$node.label[1]=100
            gt_sub$node.label[gt_sub$node.label=="0"]="1"
            gt_sub$node.label[gt_sub$node.label==""]="1"
            gt_sub$node.label[gt_sub$node.label=="NA"]="1"
            gt_sub$node.label[gt_sub$node.label=="Root"]=100
            gt_sub$node.label=as.numeric(gt_sub$node.label)/100
            rooted_phy=c(rooted_phy,write.tree(gt_sub))
               
        } 
        
    } 
    phy=rooted_phy
    for (i in 1:n_retic)
    {    
        filename=deparse(substitute(sp_list))
        f_n=paste(filename,"_phylonet_genes_",i,"ret.nex",sep="")
        write("#NEXUS\n\nBEGIN TREES;",f_n)
        write(paste("Tree fixtr = ",write.tree(sp_phy_sub)),f_n,append=TRUE)
        d=data.frame(rep("Tree",length(phy)),paste("gt",1:length(phy),"=",sep=""),phy)
        write.table(d,f_n, quote = FALSE,row.names = FALSE, col.names=FALSE,append=TRUE)
        write("END;\n\nBEGIN PHYLONET;",f_n,append=TRUE)
        write(paste("InferNetwork_ML (all)",i,"-s fixtr -di -pl 10 -x 15 -b 0.9 -po;","\nEND;"),f_n,append=TRUE)
    }    
}    


epio=c('Agriomorpha_fusca','Devadatta_kompieri','Diphlebia_euphoeoides','Euphaea_masoni','Euphaea_ochracea','Euphaea_decorata','Heteragrion_erythrogastrum','Heteragrion_majus','Miocora_notoxantha','Philogenia_carrillica','Austroargiolestes_christine','Heliocypha_perforata','Platycypha_caligata','Calopteryx_splendens','Atrocalopteryx_coomani','Mnais_costalis','Hetaerina_americana','Philoganga_vetusta','Rhinagrion_viridatum','Argia_fumipennis','Chromagrion_conditum','Nehalennia_gracilis','Mecistogaster_modesta','Megaloprepus_caerulatus','Telebasis_salva','Protoneura_sulfurata','Psaironeura_remissa','Coenagrion_puella','Cyanallagma_interruptum','Enallagma_sp','Ischnura_asiatica','Ischnura_cervula','Ischnura_verticalis','Ischnura_hastata','Ischnura_elegans','Ischnura_heterosticta','Ischnura_ramburii','Megalagrion_hawaiiense','Argiocnemis_sp','Coeliccia_sp','Copera_marginipes','Prodasineura_autumnalis','Protosticta_beaumonti','Archilestes_grandis','Indolestes_peregrinus','Episynlestes_cristatus','Synlestes_weyersii','Perissolestes_remotus','Epiophlebia_superstes','Aeshna_palmata','Anax_junius','Anax_walsinghami','Anax_parthenope','Gynacantha_tibiata','Austroaeschna_subapicalis','Telephlebia_godeffroyi','Asiagomphus_melaenops','Phanogomphus_spicatus','Stylurus_spiniceps','Leptogomphus_perforatus','Ictinogomphus_pertinax','Phenes_raptor','Tanypteryx_pryeri','Anotogaster_sieboldii','Cordulegaster_boltonii','Cordulegaster_dorsalis','Cordulegaster_maculata','Neopetalia_punctata','Chlorogomphus_auratus','Gomphomacromia_paradoxa','Eusynthemis_nigra','Macromia_amphigena','Neurocordulia_yamaskanensis','Somatochlora_uchidai','Pantala_flavescens','Rhyothemis_variegata','Ladona_fulva','Orthetrum_albistylum','Libellula_forensis','Libellula_saturata','Sympetrum_frequens','Acisoma_variegatum','Erythrodiplax_connata','Ephemera_danica','Isonychia_kiangsinensis')
anax=c('Aeshna_palmata','Anax_junius','Anax_walsinghami','Anax_parthenope','Gynacantha_tibiata','Austroaeschna_subapicalis','Telephlebia_godeffroyi',"Isonychia_kiangsinensis","Ephemera_danica") 
gompeta=c('Asiagomphus_melaenops','Phanogomphus_spicatus','Stylurus_spiniceps','Leptogomphus_perforatus','Ictinogomphus_pertinax','Phenes_raptor','Tanypteryx_pryeri','Gomphomacromia_paradoxa',"Isonychia_kiangsinensis","Ephemera_danica")
libs=c("Pantala_flavescens","Sympetrum_frequens","Rhyothemis_variegata","Erythrodiplax_connata", "Acisoma_variegatum", "Libellula_saturata","Libellula_forensis","Ladona_fulva", "Orthetrum_albistylum","Isonychia_kiangsinensis","Ephemera_danica")

phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",anax,1)
phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",gompeta,1)
phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",libs,1)
phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",epio,1) # change to MPL "InferNetwork_MPL (all) 1 -s fixtr -di -pl 15 -x 1 -h {Epiophlebia_superstes} -n 10;"