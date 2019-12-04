library(ape)
library(phytools)
library(phangorn)

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

Perilestidae="Perrisolestes_remotus"

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
fam_list=list(Calopterigoidea,Coenagrionidae,Platycnemididae,Lestidae,Synlestidae,Perilestidae,Aeshnidae,Gomphidae,Petaluridae,Cordulegastridae,Neopetaliidae,Chlorogomphidae,Synthemistidae,Macromiidae,Corduliidae,Libellulidae,Outgroup,Epiophlebiidae)
names(fam_list)=c("Calopterigoidea","Coenagrionidae","Platycnemididae","Lestidae","Synlestidae","Perilestidae","Aeshnidae","Gomphidae","Petaluridae","Cordulegastridae","Neopetaliidae","Chlorogomphidae","Synthemistidae","Macromiidae","Corduliidae","Libellulidae","Outgroup","Epiophlebiidae")
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
   
}  
all=read.tree("phylonet_trees")
d=data.frame(rep("Tree",length(all)),paste("gt",1:length(all),"=",sep=""),write.tree(all))
write.table(d,"phylonet_trees_nexus",quote = F,row.names = F, col.names=F)



fam_list=list(RZ=c(Calopterigoidea,Coenagrionidae,Platycnemididae,Platystictidae),Lestoidea=c(Lestidae,Synlestidae,Perilestidae),Aeshnidae=Aeshnidae,RA=c(Gomphidae,Petaluridae,Cordulegastridae,Neopetaliidae,Chlorogomphidae,Synthemistidae,Macromiidae,Corduliidae,Libellulidae),Outgroup=Outgroup,Epiophlebiidae=Epiophlebiidae)

fam_list=list(RZ=c(Calopterigoidea,Coenagrionidae,Platycnemididae,Platystictidae),Lestoidea=c(Lestidae,Synlestidae,Perilestidae),Aeshnidae=Aeshnidae,RA=c(Gomphidae,Petaluridae,Cordulegastridae,Neopetaliidae,Chlorogomphidae,Synthemistidae,Macromiidae,Corduliidae,Libellulidae),Epiophlebiidae=Epiophlebiidae)


get_phylonet=function(fam_list,replics)
{
    zz=1
    gene_n=0
    all_trees=c()
    full_trees=c()
    for (tr in phy)
    {
        print(zz)
        zz=zz+1
        if ("Epiophlebia_superstes" %in% tr$tip.label)
        {
            for (reps in 1:replics)
            {
                sp_intersect=lapply(fam_list,sample,1)
                tree_to_newick=keep.tip(tr,intersect(sp_intersect,tr$tip.label))
                for ( i in tree_to_newick$tip.label)
                {
                    repl=names(which(lapply(fam_list,'%in%',x=i)==T))
                    tree_to_newick$tip.label[which(i==tree_to_newick$tip.label)]=repl
                }
                all_trees=c(all_trees,write.tree(tree_to_newick))
                if (length(tree_to_newick$tip.label)==5)
                {
                    tree_to_newick=bind.tip(tree_to_newick, "Outgroup")
                    tree_to_newick=root(tree_to_newick,"Outgroup",resolve.root = T)
                    full_trees=c(full_trees,write.tree(tree_to_newick))
                }    
            }    


        }
    }    
    write.table(all_trees,"phylonet_trees",quote = F,row.names = F, col.names=F)
    write.table(full_trees,"densitree_trees",quote = F,row.names = F, col.names=F)
    all=read.tree("phylonet_trees")
    write.table("#NEXUS\n\nBEGIN TREES;","phylonet_trees_nexus",quote = F,row.names = F, col.names=F)
    write.table("Tree fixtr = ((Lestoidea,RZ),(Epiophlebiidae,(Aeshnidae,RA)));","phylonet_trees_nexus",quote = F,row.names = F, col.names=F,append=T)
    d=data.frame(rep("Tree",length(all)),paste("gt",1:length(all),"=",sep=""),write.tree(all))
    write.table(d,"phylonet_trees_nexus",quote = F,row.names = F, col.names=F,append=T)
    write.table("END;\n\nBEGIN PHYLONET;","phylonet_trees_nexus",quote = F,row.names = F, col.names=F,append=T)
    l=length(all)
    seq_left=seq(1,l,by=replics)
    seq_right=seq(0,l,by=replics)
    pp=data.frame(b=paste("{gt",seq_left,sep=""),a=paste("-gt",seq_right[2:length(seq_right)],"}",sep=""))
    parts=paste(apply(pp,1,paste,collapse=""),collapse=",")
    write.table(paste("InferNetwork_MPL","(",parts,")",1,"-s fixtr -fs -di -pl 15;","\nEND;"),"phylonet_trees_nexus",quote = F,row.names = F, col.names=F,append=T)   
}





###All trees
zz=read.tree("/Users/Anton/Downloads/trees/all.trees")
fordensitree=function(alltree)
{
    all_trees=c()

    for (tr in alltree)
    {
        tr=root(tr,"Isonychia_kiangsinensis",resolve.root = T)
        if ("Ephemera_danica" %in% tr$tip.label)
        {
            tr=drop.tip(tr,"Ephemera_danica")
        }
        tr=drop.tip(tr,"Isonychia_kiangsinensis")
        tr$node.label=NULL
        tr$edge.length=NULL
        all_trees=c(all_trees,write.tree(tr))
        
    }
    write.table(all_trees,"fullphylogeny_all",quote = F,row.names = F, col.names=F)
}    



### All trees summary 
kk=read.tree("/Users/Anton/Downloads/fullphylogeny_all")



sporder=c("Perissolestes_remotus","Synlestes_weyersii","Episynlestes_cristatus","Indolestes_peregrinus","Archilestes_grandis","Protosticta_beaumonti","Euphaea_decorata","Euphaea_ochracea","Euphaea_masoni","Diphlebia_euphoeoides","Devadatta_kompieri","Agriomorpha_fusca","Philogenia_carrillica","Miocora_notoxantha","Heteragrion_majus","Heteragrion_erythrogastrum","Hetaerina_americana","Mnais_costalis","Atrocalopteryx_coomani","Calopteryx_splendens","Platycypha_caligata","Heliocypha_perforata","Austroargiolestes_christine","Rhinagrion_viridatum","Philoganga_vetusta","Prodasineura_autumnalis","Copera_marginipes","Coeliccia_sp","Telebasis_salva","Megaloprepus_caerulatus","Mecistogaster_modesta","Nehalennia_gracilis","Chromagrion_conditum","Psaironeura_remissa","Protoneura_sulfurata","Argia_fumipennis","Coenagrion_puella","Argiocnemis_sp","Megalagrion_hawaiiense","Ischnura_ramburii","Ischnura_heterosticta","Ischnura_elegans","Ischnura_hastata","Ischnura_verticalis","Ischnura_cervula","Ischnura_asiatica","Enallagma_sp","Cyanallagma_interruptum","Epiophlebia_superstes","Telephlebia_godeffroyi","Austroaeschna_subapicalis","Gynacantha_tibiata","Anax_parthenope","Anax_walsinghami","Anax_junius","Aeshna_palmata","Tanypteryx_pryeri","Phenes_raptor","Ictinogomphus_pertinax","Leptogomphus_perforatus","Stylurus_spiniceps","Phanogomphus_spicatus","Asiagomphus_melaenops","Chlorogomphus_auratus","Neopetalia_punctata","Cordulegaster_maculata","Cordulegaster_dorsalis","Cordulegaster_boltonii","Anotogaster_sieboldii","Eusynthemis_nigra","Gomphomacromia_paradoxa","Macromia_amphigena","Somatochlora_uchidai","Neurocordulia_yamaskanensis","Rhyothemis_variegata","Pantala_flavescens","Libellula_saturata","Libellula_forensis","Orthetrum_albistylum","Ladona_fulva","Sympetrum_frequens","Erythrodiplax_connata","Acisoma_variegatum")


quartz(width=8.21, height=10)
densiTree(kk,alpha = 0.2,scaleX = T,jitter = list(amount = 0.1, random=TRUE),col=rep(c("blue","red","blue","red","blue","darkgreen","blue","red","blue","red","blue","red"),c(13,4,3,1,1,1,1,2,1,4,13,4)),consensus=rev(sporder),tip.color="black",label.offset=0.01,cex=0.6,scale.bar = F)
legend("topleft",c("DNA","Protein","AF") ,lty=c(1,1,1),col=c("blue","red","darkgreen"),text.col ="black",lwd=3)
quartz.save("All_supermatrix_trees.jpeg", type = "jpeg",antialias=F,bg="white",dpi=400,pointsize=12)


#Gene trees densitree
genefordensitree=function(alltree,filename)
{
    all_trees=c()

    for (tr in alltree)
    {
        if (any(c("Ephemera_danica","Isonychia_kiangsinensis") %in% tr$tip.label))
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
        
            tr$node.label=NULL
            tr$edge.length=NULL
            all_trees=c(all_trees,write.tree(tr))
        }
    }
    write.table(all_trees,filename,quote = F,row.names = F, col.names=F)
}    

filter_by_n=function(alltree,g1,g2,g3,ng1,ng2,ng3,filename)
{
    all_trees=c()
    for (tr in alltree)
    {
        l1=length(intersect(tr$tip.label,g1))
        l2=length(intersect(tr$tip.label,g2))
        l3=length(intersect(tr$tip.label,g3))
        if (l1 >= ng1 & l2 >= ng2 & l3 >= ng3)
        {
           all_trees=c(all_trees,write.tree(tr)) 
        }    
    } 
    write.table(all_trees,filename,quote = F,row.names = F, col.names=F)
}    

filter_by_boot=function(alltree,bootcut,filename)
{
    all_trees=c()
    for (tr in alltree)
    {
        meanboot=mean(as.numeric(tr$node.label[tr$node.label!=""]))
        if (meanboot >= bootcut)
        {
           all_trees=c(all_trees,write.tree(tr)) 
        }    
    } 
    write.table(all_trees,filename,quote = F,row.names = F, col.names=F)
}    

g1=c("Perissolestes_remotus","Synlestes_weyersii","Episynlestes_cristatus","Indolestes_peregrinus","Archilestes_grandis","Protosticta_beaumonti","Euphaea_decorata","Euphaea_ochracea","Euphaea_masoni","Diphlebia_euphoeoides","Devadatta_kompieri","Agriomorpha_fusca","Philogenia_carrillica","Miocora_notoxantha","Heteragrion_majus","Heteragrion_erythrogastrum","Hetaerina_americana","Mnais_costalis","Atrocalopteryx_coomani","Calopteryx_splendens","Platycypha_caligata","Heliocypha_perforata","Austroargiolestes_christine","Rhinagrion_viridatum","Philoganga_vetusta","Prodasineura_autumnalis","Copera_marginipes","Coeliccia_sp","Telebasis_salva","Megaloprepus_caerulatus","Mecistogaster_modesta","Nehalennia_gracilis","Chromagrion_conditum","Psaironeura_remissa","Protoneura_sulfurata","Argia_fumipennis","Coenagrion_puella","Argiocnemis_sp","Megalagrion_hawaiiense","Ischnura_ramburii","Ischnura_heterosticta","Ischnura_elegans","Ischnura_hastata","Ischnura_verticalis","Ischnura_cervula","Ischnura_asiatica","Enallagma_sp","Cyanallagma_interruptum")

g2="Epiophlebia_superstes"

g3=c("Telephlebia_godeffroyi","Austroaeschna_subapicalis","Gynacantha_tibiata","Anax_parthenope","Anax_walsinghami","Anax_junius","Aeshna_palmata","Tanypteryx_pryeri","Phenes_raptor","Ictinogomphus_pertinax","Leptogomphus_perforatus","Stylurus_spiniceps","Phanogomphus_spicatus","Asiagomphus_melaenops","Chlorogomphus_auratus","Neopetalia_punctata","Cordulegaster_maculata","Cordulegaster_dorsalis","Cordulegaster_boltonii","Anotogaster_sieboldii","Eusynthemis_nigra","Gomphomacromia_paradoxa","Macromia_amphigena","Somatochlora_uchidai","Neurocordulia_yamaskanensis","Rhyothemis_variegata","Pantala_flavescens","Libellula_saturata","Libellula_forensis","Orthetrum_albistylum","Ladona_fulva","Sympetrum_frequens","Erythrodiplax_connata","Acisoma_variegatum")



buscoprot=read.tree("/Users/Anton/Downloads/gene_trees/BUSCO50_prot_pasta_iqtree_all")
genefordensitree(buscoprot,"buscoprot_densitree")
buscoprot_tr=read.tree("/Users/Anton/Downloads/buscoprot_densitree")

buscodna=read.tree("/Users/Anton/Downloads/gene_trees/BUSCO50_dna_pasta_iqtree_all")
genefordensitree(buscodna,"buscodna_densitree")
buscodna_tr=read.tree("/Users/Anton/Downloads/buscodna_densitree")
densiTree(c(kk[3],buscodna_tr),scaleX = T,consensus=rev(sporder),jitter = list(amount = 0.2, random=TRUE),col=c(adjustcolor("navy", alpha.f = 0.2),adjustcolor( rep("black",length(buscodna_tr)), alpha.f = 0.01)),alpha=1,tip.color="black",label.offset=0.01,cex=0.6,scale.bar = F)
legend("topleft",title="BUSCO50",c("species tree","gene tree") ,lty=c(1,1),col=c("red","black"),text.col ="black",lwd=3)

yangdna=read.tree("/Users/Anton/Downloads/gene_trees/YANG50_dna_pasta_iqtree_all")
genefordensitree(yangdna,"yangdna_densitree")
yangdna_tr=read.tree("/Users/Anton/Downloads/yangdna_densitree")
densiTree(c(kk[3],yangdna_tr),scaleX = T,consensus=rev(sporder),jitter = list(amount = 0.2, random=TRUE),col=c(adjustcolor("navy", alpha.f = 0.5),adjustcolor( rep("black",length(yangdna_tr)), alpha.f = 0.007)),alpha=1,tip.color="black",label.offset=0.01,cex=0.6,scale.bar = F)
legend("topleft",title="YANG50",c("species tree","gene tree") ,lty=c(1,1),col=c("red","black"),text.col ="black",lwd=3)

yangprot=read.tree("/Users/Anton/Downloads/gene_trees/YANG50_prot_pasta_iqtree_all")
genefordensitree(yangprot,"yangprot_densitree")
yangprot_tr=read.tree("/Users/Anton/Downloads/yangprot_densitree")


densiTree(buscoprot_tr,alpha = 0.05,scaleX = T,consensus=rev(sporder))



par(bg = 'black')
densiTree(c(yangdna_tr,kk[3]),scaleX = T,consensus=rev(sporder),jitter = list(amount = 0.3, random=TRUE),col=c(adjustcolor( rep("white",length(yangdna_tr)), alpha.f = 0.007),"red"),alpha=1,lwd=4)