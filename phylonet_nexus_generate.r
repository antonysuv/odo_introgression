library('ape')
library('phytools')
library('phangorn')



#Listselect
listsel=function(l,ind)
{
    new_l=l
    for (i in 1:length(l))
    new_l[[i]]=l[[i]][ind[[i]]]
    return(new_l)
        
}    

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
        if (any(sp_list %in% gt$tip.label) & any(out_v %in% gt$tip.label) & sum(sp_list %in% gt$tip.label)>=3)
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
        write(paste("InferNetwork_MPL (all)",i,"-s fixtr -di -pl 10 -x 100 -b 0.9 -n 5 -po;","\nEND;"),f_n,append=TRUE)
        #return(phy)
    }    
}    


epio=c('Ischnura_elegans','Copera_marginipes','Protosticta_beaumonti','Archilestes_grandis','Indolestes_peregrinus','Episynlestes_cristatus','Synlestes_weyersii','Perissolestes_remotus','Epiophlebia_superstes','Aeshna_palmata','Anax_junius','Anax_walsinghami','Anax_parthenope','Gynacantha_tibiata','Austroaeschna_subapicalis','Telephlebia_godeffroyi','Phenes_raptor','Ladona_fulva','Ephemera_danica','Isonychia_kiangsinensis')
anax=c('Aeshna_palmata','Anax_junius','Anax_walsinghami','Anax_parthenope','Gynacantha_tibiata','Austroaeschna_subapicalis','Telephlebia_godeffroyi',"Isonychia_kiangsinensis","Ephemera_danica") 
gompeta=c("Anax_parthenope",'Asiagomphus_melaenops','Phanogomphus_spicatus','Stylurus_spiniceps','Leptogomphus_perforatus','Ictinogomphus_pertinax','Phenes_raptor','Tanypteryx_pryeri','Gomphomacromia_paradoxa',"Isonychia_kiangsinensis","Ephemera_danica")
libs=c("Pantala_flavescens","Sympetrum_frequens","Rhyothemis_variegata","Erythrodiplax_connata", "Acisoma_variegatum", "Libellula_saturata","Libellula_forensis","Ladona_fulva", "Orthetrum_albistylum","Isonychia_kiangsinensis","Ephemera_danica")

phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",anax,1)
phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",gompeta,1)
phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",libs,1)
phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",epio,1) # change to MPL "InferNetwork_MPL (all) 1 -s fixtr -di -pl 15 -x 1 -h {Epiophlebia_superstes} -n 10;"



SuperMatrix_50BUSCO_dna_pasta_ali_trim





#ML Phylonet species random sample
phylonet_chunk_topo_random=function(sp_tree,gene_trees,fam_list,file_n,reps)
{
    sp_phy=read.tree(sp_tree)
    sp_phy$node.label=NULL
    out_v=c("Isonychia_kiangsinensis","Ephemera_danica")
    sp_list=lapply(fam_list,sample,1)
    sp_phy_sub=keep.tip(sp_phy,as.character(sp_list[sp_list %in% sp_phy$tip.label]))
    for ( i in sp_phy_sub$tip.label)
    {
        repl=names(which(lapply(fam_list,'%in%',x=i)==T))
        sp_phy_sub$tip.label[which(i==sp_phy_sub$tip.label)]=repl
    }
    phy=read.tree(gene_trees)
    rooted_phy=c()
    for (gt in phy)
    {
       
       if (all(unlist(lapply(lapply(fam_list,"%in%",gt$tip.label),any))))
       {    
            for (re in 1:reps)
            {
                pos=lapply(lapply(lapply(fam_list,"%in%",gt$tip.label),which),sample,1)
                sp_list=listsel(fam_list,pos)
                if (all(sp_list %in% gt$tip.label) & any(out_v %in% gt$tip.label))
                {

                    out_sp=out_v[which(out_v %in% gt$tip.label)[1]]
                    gt_r=root(gt,out_sp,resolve.root=T)
                    gt_sub=keep.tip(gt_r,unlist(sp_list))
                    gt_sub$node.label[1]=100
                    gt_sub$node.label[gt_sub$node.label=="0"]="1"
                    gt_sub$node.label[gt_sub$node.label==""]="1"
                    gt_sub$node.label[gt_sub$node.label=="NA"]="1"
                    gt_sub$node.label[gt_sub$node.label=="Root"]=100
                    gt_sub$node.label=as.numeric(gt_sub$node.label)/100
                    for ( i in gt_sub$tip.label)
                    {
                        repl=names(which(lapply(fam_list,'%in%',x=i)==T))
                        gt_sub$tip.label[which(i==gt_sub$tip.label)]=repl
                    }
                    rooted_phy=c(rooted_phy,write.tree(gt_sub))

                } 

            }      
        }
    } 
    phy=rooted_phy  
    filename=deparse(substitute(fam_list))
    f_n=paste(filename,"_phylonet_genes_",file_n,"MCMC.nex",sep="")
    write("#NEXUS\n\nBEGIN TREES;",f_n)
    write(paste("Tree fixtr = ",write.tree(sp_phy_sub)),f_n,append=TRUE)
    d=data.frame(rep("Tree",length(phy)),paste("gt",1:length(phy),"=",sep=""),phy)
    write.table(d,f_n, quote = FALSE,row.names = FALSE, col.names=FALSE,append=TRUE)
    write("END;\n\nBEGIN PHYLONET;",f_n,append=TRUE)
    write(paste("MCMC_GT (all) -cl 1010000 -bl 10000 -sf 1000 -mr 1 -pl 10;","\nEND;"),f_n,append=TRUE)
        #return(phy)
      
}    



#Subsample Epiophlebiidae
Lestoidea=c("Perissolestes_remotus","Synlestes_weyersii","Episynlestes_cristatus","Indolestes_peregrinus","Archilestes_grandis","Protosticta_beaumonti")
RZ=c("Euphaea_decorata","Euphaea_ochracea","Euphaea_masoni","Diphlebia_euphoeoides","Devadatta_kompieri","Agriomorpha_fusca","Philogenia_carrillica","Miocora_notoxantha","Heteragrion_majus","Heteragrion_erythrogastrum","Hetaerina_americana","Mnais_costalis","Atrocalopteryx_coomani","Calopteryx_splendens","Platycypha_caligata","Heliocypha_perforata","Austroargiolestes_christine","Rhinagrion_viridatum","Philoganga_vetusta","Prodasineura_autumnalis","Copera_marginipes","Coeliccia_sp","Telebasis_salva","Megaloprepus_caerulatus","Mecistogaster_modesta","Nehalennia_gracilis","Chromagrion_conditum","Psaironeura_remissa","Protoneura_sulfurata","Argia_fumipennis","Coenagrion_puella","Argiocnemis_sp","Megalagrion_hawaiiense","Ischnura_ramburii","Ischnura_heterosticta","Ischnura_elegans","Ischnura_hastata","Ischnura_verticalis","Ischnura_cervula","Ischnura_asiatica","Enallagma_sp","Cyanallagma_interruptum")
Epiophlebiidae="Epiophlebia_superstes"
Aeshnidae=c("Telephlebia_godeffroyi", "Austroaeschna_subapicalis","Gynacantha_tibiata","Anax_parthenope","Anax_walsinghami","Anax_junius","Aeshna_palmata")
RA=c("Tanypteryx_pryeri","Phenes_raptor","Ictinogomphus_pertinax","Leptogomphus_perforatus","Stylurus_spiniceps","Phanogomphus_spicatus","Asiagomphus_melaenops","Chlorogomphus_auratus","Neopetalia_punctata","Cordulegaster_maculata","Cordulegaster_dorsalis","Cordulegaster_boltonii","Anotogaster_sieboldii","Eusynthemis_nigra","Gomphomacromia_paradoxa","Macromia_amphigena","Somatochlora_uchidai","Neurocordulia_yamaskanensis","Rhyothemis_variegata","Pantala_flavescens","Libellula_saturata","Libellula_forensis","Orthetrum_albistylum","Ladona_fulva","Sympetrum_frequens","Erythrodiplax_connata","Acisoma_variegatum")
Outgroup=c("Isonychia_kiangsinensis","Ephemera_danica")
epio_run=list(Lestoidea=Lestoidea,Epiophlebiidae=Epiophlebiidae,Aeshnidae=Aeshnidae,RA=RA,RZ=RZ,Outgroup=Outgroup)

for (i in 1:5)
{
    phylonet_chunk_topo_random("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",epio_run,i,10)
}


#Subsample Gomphidae-Petaluridae
RA1=c("Epiophlebia_superstes","Telephlebia_godeffroyi", "Austroaeschna_subapicalis","Gynacantha_tibiata","Anax_parthenope","Anax_walsinghami","Anax_junius","Aeshna_palmata")
Petaluridae=c("Tanypteryx_pryeri","Phenes_raptor")
Gomphidae=c("Ictinogomphus_pertinax","Leptogomphus_perforatus","Stylurus_spiniceps","Phanogomphus_spicatus","Asiagomphus_melaenops")
RA2=c("Chlorogomphus_auratus","Neopetalia_punctata","Cordulegaster_maculata","Cordulegaster_dorsalis","Cordulegaster_boltonii","Anotogaster_sieboldii","Eusynthemis_nigra","Gomphomacromia_paradoxa","Macromia_amphigena","Somatochlora_uchidai","Neurocordulia_yamaskanensis","Rhyothemis_variegata","Pantala_flavescens","Libellula_saturata","Libellula_forensis","Orthetrum_albistylum","Ladona_fulva","Sympetrum_frequens","Erythrodiplax_connata","Acisoma_variegatum")
Outgroup=c("Isonychia_kiangsinensis","Ephemera_danica")
gompeta_run=list(Petaluridae=Petaluridae,Gomphidae=Gomphidae,RA1=RA1,RA2=RA2,Outgroup=Outgroup)

#MCMC
for (i in 1:5)
{
    phylonet_chunk_topo_random("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",gompeta_run,i,10)
}    



#ML
phylonet_chunk_topo_random("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",gompeta_run,1,10)



log_epio=c()
for (i in mm)
{
    ww="Epiophlebia_superstes" %in% i$tip.label
    log_epio=c(log_epio,ww)
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
#Epiophlebia test 
















anax_sub=c('Aeshna_palmata','Anax_junius','Anax_walsinghami','Anax_parthenope','Gynacantha_tibiata','Austroaeschna_subapicalis','Telephlebia_godeffroyi')
anax_tt=read.tree(text=phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",anax_sub,1))

densiTree(c(anax_tt,anax_tt[3]),scaleX = F,consensus=anax_sub,col=c(adjustcolor("navy", alpha.f = 0.2),adjustcolor( rep("black",length(anax_tt)), alpha.f = 0.01)),alpha=1,tip.color="black",label.offset=0.01,cex=0.6,scale.bar = F)

gompeta_sub=c('Asiagomphus_melaenops','Phanogomphus_spicatus','Stylurus_spiniceps','Leptogomphus_perforatus','Ictinogomphus_pertinax','Phenes_raptor','Tanypteryx_pryeri','Gomphomacromia_paradoxa')
gompeta_tt=read.tree(text=phylonet_chunk_topo("BUSCO50_dna_pasta_nopart_iqtree_root.tre","BUSCO50_dna_pasta_iqtree_all_wboot",gompeta_sub,1))
densiTree(gompeta_tt,scaleX = F,consensus=anax_sub,col=c(adjustcolor("navy", alpha.f = 0.2),adjustcolor( rep("black",length(anax_tt)), alpha.f = 0.01)),alpha=1,tip.color="black",label.offset=0.01,cex=0.6,scale.bar = F)



