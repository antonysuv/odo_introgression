library(ape)

Calopterigoidea=c("Agriomorpha_fusca",
	"Devedatta_sp",
	"Diphlabia_sp",
	"Euphaea_masoni",
	"Euphaea_sp",
	"Euphaea_sp2",
	"Heteragrion_erythrogastrum",
	"Heteragrion_majus",
	"Polythore_notoxantha",
	"Philogenia_carillaca",
	"Austroargiolestes_christine",
	"Libellago_sp",
	"Platychypha_sp",
	"Calopteryx_splendens",
	"Vestalis_sp",
	"Mnais_costalis",
	"Hetaerina_americana_milly",
	"Philoganga_vetusta",
	"Rhinoagrion_sp")

Coenagrionidae=c("Argia_fumipennis",
	"Chromagrion_conditum",
	"Nehalennia_gracilis",
	"Mecistogaster_modesta",
	"Megaloprepus_caerulatus",
	"Telebasis_salva",
	"Protoneura_sulfurata",
	"Psaironeura_remissa",
	"Coenagrion_puella",
	"Cyanallagma_interruptum",
	"Enallagma_sp",
	"Ischnura_asiatica",
	"Ischnura_cervula",
	"Ischnura_verticalis",
	"Ischnura_hastata",
	"Ischnura_elegans",
	"Ischnura_heterosticta",
	"Ischnura_ramburii",
	"Megalagrion_hawaiiensis",
	"Ellatoneura_sp")

Platycnemididae=c("Coeliccia_sp",
	"Copera_manginipes",
	"Prodasineura_automalis")
	
Platystictidae="Protosticta_beaumonti"

Lestidae=c("Archilestes_grandis",
	"Indolestes_peregrinus")
	
Synlestidae=c("Episylestes_cristatus",
	"Synlestes_weyersii")

Preilestidae="Perrisolestes_remotus"

Epiophlebiidae="Epiophlebia_superstes"

Aeshnidae=c("Aeshna_palmata",
	"Anax_junius",
	"Anax_walsinghamii",
	"Anax_parthenope",
	"Gynacantha_tibiata",
	"Austroaeshna_subapicalis",
	"Telephlebia_sp")

Gomphidae=c("Asiagomphus_melaenops",
	"Gomphus_spicatus",
	"Stylurus_spiniceps",
	"Leptogomphus_perforatus",
	"Ictinogomphus_sp")

Petaluridae=c("Phenes_raptor",
	"Tanypteryx_pryerii")
	
Cordulegastridae=c("Anotogaster_sieboldii",
	"Cordulegaster_boltonii",
	"Cordulegaster_dorsalis",
	"Cordulegaster_maculata")

Neopetaliidae="Neopetalia_punctata"

Chlorogomphidae="Chlorogomphus_auratus"

Synthemistidae=c("Gomphomarcromia_paradoxa",
	"Synthemis_sp")
Macromiidae="Macromia_amphigena"

Corduliidae=c("Neurocordulia_yamaskanensis",
	"Somatochlora_uchidai")

Epiophlebiidae="Epiophlebia_superstes"

Libellulidae=c("Pantala_flavscens",
	"Rhyothemis_sp",
	"Ladona_fulva",
	"Orthetrum_albistylum",
	"Libellula_forensis",
	"Libellula_saturnata",
	"Sympetrum_frequens",
	"Acisoma_variegatum",
	"Erythrodiplax_conata")
Outgroup=c("Ephemera_danica",
	"Isonychia_kiangsinensis")

phy=read.tree("/Users/Anton/Downloads/50BUSCO_dna_pasta_iqtree_trees_all")
fam_list=list(Calopterigoidea,Coenagrionidae,Platycnemididae,Lestidae,Synlestidae,Preilestidae,Aeshnidae,Gomphidae,Petaluridae,Cordulegastridae,Neopetaliidae,Chlorogomphidae,Synthemistidae,Macromiidae,Corduliidae,Libellulidae,Outgroup,Epiophlebiidae)
names(fam_list)=c("Calopterigoidea","Coenagrionidae","Platycnemididae","Lestidae","Synlestidae","Preilestidae","Aeshnidae","Gomphidae","Petaluridae","Cordulegastridae","Neopetaliidae","Chlorogomphidae","Synthemistidae","Macromiidae","Corduliidae","Libellulidae","Outgroup","Epiophlebiidae")
for (tr in phy)
{
    tip_list=c()
    if (any(tr$tip.label %in% "Epiophlebia_superstes"))
    {
        for (pos_fam in 1:length(fam_list))
        {
            fam=unlist(fam_list[pos_fam])
            if(any(tr$tip.label %in% fam))
            {
                select=sample(tr$tip.label[tr$tip.label %in% fam],1)
                tip_list=c(tip_list,select)
            }    
        }       
        tree_to_newick=keep.tip(tr,tip_list)
        for ( i in tree_to_newick$tip.label)
        {
            repl=names(which(lapply(fam_list,'%in%',x=i)==T))
            tree_to_newick$tip.label[which(i==tree_to_newick$tip.label)]=repl
        }    
        write.tree(tree_to_newick,"phylonet_trees",append=T)

    }
    all=read.tree("phylonet_trees")
    d=data.frame(rep("Tree",length(all)),paste("gt",1:length(all),"=",sep=""),write.tree(all))
    write.table(d,"phylonet_trees_nexus",quote = F,row.names = F, col.names=F)
}    