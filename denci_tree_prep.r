library(ape)
library(phytools)
library(phangorn)

genefordensitree=function(alltree)
{
    all_trees=c()

    for (tr in alltree)
    {
        if (any(c("Ephemera_danica","Isonychia_kiangsinensis","Outgroup") %in% tr$tip.label) & length(tr$tip.label) > 3)
        {
            if ("Ephemera_danica" %in% tr$tip.label)
            {
                tr=root(tr,"Ephemera_danica",resolve.root = T)
                tr=drop.tip(tr,"Ephemera_danica")
            }  
            
            if ("Isonychia_kiangsinensis" %in% tr$tip.label)
            {
                tr=root(tr,"Isonychia_kiangsinensis",resolve.root = T)
                tr=drop.tip(tr,"Isonychia_kiangsinensis")
            } 
            
            if ("Outgroup" %in% tr$tip.label)
            {
                tr=root(tr,"Outgroup",resolve.root = T)
                #tr=drop.tip(tr,"Outgroup")
            } 
        
            tr$node.label=NULL
            tr$edge.length=NULL
            all_trees=c(all_trees,write.tree(tr))
        }
    }
    return(read.tree(text=all_trees))
}    


sporder=c("Perissolestes_remotus","Synlestes_weyersii","Episynlestes_cristatus","Indolestes_peregrinus","Archilestes_grandis","Protosticta_beaumonti","Euphaea_decorata","Euphaea_ochracea","Euphaea_masoni","Diphlebia_euphoeoides","Devadatta_kompieri","Agriomorpha_fusca","Philogenia_carrillica","Miocora_notoxantha","Heteragrion_majus","Heteragrion_erythrogastrum","Hetaerina_americana","Mnais_costalis","Atrocalopteryx_coomani","Calopteryx_splendens","Platycypha_caligata","Heliocypha_perforata","Austroargiolestes_christine","Rhinagrion_viridatum","Philoganga_vetusta","Prodasineura_autumnalis","Copera_marginipes","Coeliccia_sp","Telebasis_salva","Megaloprepus_caerulatus","Mecistogaster_modesta","Nehalennia_gracilis","Chromagrion_conditum","Psaironeura_remissa","Protoneura_sulfurata","Argia_fumipennis","Coenagrion_puella","Argiocnemis_sp","Megalagrion_hawaiiense","Ischnura_ramburii","Ischnura_heterosticta","Ischnura_elegans","Ischnura_hastata","Ischnura_verticalis","Ischnura_cervula","Ischnura_asiatica","Enallagma_sp","Cyanallagma_interruptum","Epiophlebia_superstes","Telephlebia_godeffroyi","Austroaeschna_subapicalis","Gynacantha_tibiata","Anax_parthenope","Anax_walsinghami","Anax_junius","Aeshna_palmata","Tanypteryx_pryeri","Phenes_raptor","Ictinogomphus_pertinax","Leptogomphus_perforatus","Stylurus_spiniceps","Phanogomphus_spicatus","Asiagomphus_melaenops","Chlorogomphus_auratus","Neopetalia_punctata","Cordulegaster_maculata","Cordulegaster_dorsalis","Cordulegaster_boltonii","Anotogaster_sieboldii","Eusynthemis_nigra","Gomphomacromia_paradoxa","Macromia_amphigena","Somatochlora_uchidai","Neurocordulia_yamaskanensis","Rhyothemis_variegata","Pantala_flavescens","Libellula_saturata","Libellula_forensis","Orthetrum_albistylum","Ladona_fulva","Sympetrum_frequens","Erythrodiplax_connata","Acisoma_variegatum")


### All species trees summary 
all_sptrees=read.tree("/Users/Anton/Downloads/fullphylogeny_all")
quartz(width=8.21, height=10)
densiTree(all_sptrees,alpha = 0.2,scaleX = T,jitter = list(amount = 0.1, random=TRUE),col=rep(c("blue","red","blue","red","blue","darkgreen","blue","red","blue","red","blue","red"),c(13,4,3,1,1,1,1,2,1,4,13,4)),consensus=rev(sporder),tip.color="black",label.offset=0.01,cex=0.6,scale.bar = F)


### All gene trees summary 
all_genes=read.tree("/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all")
quartz(width=8.21, height=10)
all_genes_trim=genefordensitree(all_genes) #1475 phylogenetic trees
densiTree(all_genes_trim,alpha = 0.009,scaleX = T,jitter = list(amount = 0.1, random=TRUE),consensus=rev(sporder),col="black",label.offset=0.01,cex=0.6,scale.bar = F)

