library("ape")
library('ggplot2')
library('pals')
library('reshape2')
library('dplyr')
library('gridExtra')


.write.tree2 <- function(phy, digits = 10, tree.prefix = "", check_tips)
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    #if (check_tips) phy$tip.label <- checkLabel(phy$tip.label)
    #if (nodelab) phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, "f", sep = "")
    cp <- function(x){
        STRING[k] <<- x
        k <<- k + 1
    }
    add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
            if (j > n) add.internal(j)
            else add.terminal(ind[j])
            if (j != desc[length(desc)]) cp(",")
        }
        cp(")")
        if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[ind[i]]))
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i]))
        }
    }

    n <- length(phy$tip.label)

    ## borrowed from phangorn:
    parent <- phy$edge[, 1]
    children <- phy$edge[, 2]
    kids <- vector("list", n + phy$Nnode)
    for (i in 1:length(parent))
        kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])

    ind <- match(1:max(phy$edge), phy$edge[, 2])

    LS <- 4*n + 5
    if (brl) LS <- LS + 4*n
    if (nodelab)  LS <- LS + n
    STRING <- character(LS)
    k <- 1
    cp(tree.prefix)
    cp("(")
    getRoot <- function(phy)
        phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
    root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
    desc <- kids[[root]]
    for (j in desc) {
        if (j > n) add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) cp(",")
    }

    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    paste(STRING, collapse = "")
}

assignInNamespace(".write.tree2", .write.tree2, "ape")



quible_trees=function(triplet,gene_trees,file_id)
{
    phy=read.tree(gene_trees)
    trees=c()
    out_v=c("Isonychia_kiangsinensis","Ephemera_danica")
    for (gt in phy)
    {
        if (all(triplet %in% gt$tip.label) & any(out_v %in% gt$tip.label))
        {
            out_sp=out_v[which(out_v %in% gt$tip.label)[1]]
            gt$tip.label[gt$tip.label==out_sp]="Outgroup"
            gt$node.label=NULL
            triplet_keep=gt$tip.label[gt$tip.label %in% triplet]
            gt_sub=keep.tip(gt,c(triplet_keep,"Outgroup"))
            trees=c(trees,write.tree(gt_sub))
        }    
        
    }    
    config="
[Input]
treefile: ./INNAME
numdistributions: 2
likelihoodthresh: 0.01
numsteps: 50
gradascentscalar: 0.5
totaloutgroup: Outgroup
multiproc: True
maxcores:1000
[Output]
OutputPath: ./OUTNAME"    
    
    config=gsub("INNAME",paste("quibltrees_",file_id,sep=""),config)
    config=gsub("OUTNAME",paste("quibl_out",file_id,sep=""),config)
    write(config,paste("quiblin_",file_id,sep=""))
    write.table(trees,paste("quibltrees_",file_id,sep=""), quote = FALSE,row.names = FALSE, col.names=FALSE,append=TRUE)
    
}    



all_triplets=function(taxa_list,gene_trees,dir_name)
{
    dir.create(dir_name)
    setwd(dir_name)
    taxa_combn=combn(taxa_list,m=3)
    for (i in 1:ncol(taxa_combn))
    {
        triplet=taxa_combn[,i]
        quible_trees(triplet,gene_trees,i)
    }    
    setwd("..")
}    


get_intropair=function(m,target_sp)
{
    l=m[,c("P1","P2","P3")]!=m[,c("outgroup","outgroup","outgroup")]
    sp_list=m[,c("P1","P2","P3")]
    pair_v=c()
    for (i in 1:nrow(l))        
    {
        pair_intro=as.character(sp_list[i,][l[i,]])
        pair_intro=paste(pair_intro[!pair_intro %in% target_sp],target_sp)
        pair_v=c(pair_v,pair_intro)
    } 
    return(pair_v)
}








tt=read.tree("BUSCO50_dna_pasta_nopart_iqtree_root.tre")

Anisozygoptera=c(extract.clade(tt,131)$tip.label,extract.clade(tt,137)$tip.label,"Epiophlebia_superstes")
Aeshnoidea=c(extract.clade(tt,137)$tip.label,extract.clade(tt,144)$tip.label)
Calopterygoidea=extract.clade(tt,91)$tip.label
Coenagrionoidea=extract.clade(tt,109)$tip.label
Lestoidea=extract.clade(tt,131)$tip.label
Cordulegastroidea=extract.clade(tt,151)$tip.label
Libelluloidea=extract.clade(tt,156)$tip.label
Zygoptera=c(Calopterygoidea,Coenagrionoidea,Lestoidea)
Anisozygoptera=c("Epiophlebia_superstes")
Anisoptera=c(Aeshnoidea,Cardulegastroidea,Libelluloidea)

Aeshnidae=extract.clade(tt,137)$tip.label
Gomphidae_Petaluridae=extract.clade(tt,144)$tip.label
Libellulidae=extract.clade(tt,161)$tip.label


all_triplets(Aeshnoidea,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Aeshnoidea_quibl")
all_triplets(Calopterygoidea,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Calopterygoidea_quibl")
all_triplets(Coenagrionoidea,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Coenagrionoidea_quibl")
all_triplets(Lestoidea,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Lestoidea_quibl")
all_triplets(Cardulegastroidea,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Cardulegastroidea_quibl")
all_triplets(Libelluloidea,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Libelluloidea_quibl")
all_triplets(Anisozygoptera,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Anisozygoptera_quibl")



################################################################ RUN QuiBL ########################################################################


total=read.csv("/Users/Anton/Downloads/quibl_all.txt",stringsAsFactors=FALSE)
P1=gsub(" ","_",paste(unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",1)),unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",2))))
P2=gsub(" ","_",paste(unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",3)),unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",4))))
P3=gsub(" ","_",paste(unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",5)),unlist(lapply(strsplit(as.character(total$triplet), "_"),"[",6))))
total$P1=P1
total$P2=P2
total$P3=P3
total=total[!(!apply(apply(total[,c("P1","P2","P3")],2,"%in%","Epiophlebia_superstes"),1,any) & total$superfamily=="Anisozygoptera"),]
total$common=FALSE
for (trip in unique(as.character(total$triplet)))
{
  maxVal=max(subset(total,triplet==trip)$count)
  # Handle ties
  if(nrow(total[which(total$triplet==trip & total$count==maxVal),])>1)
  {
     total[which(total$triplet==trip & total$count==maxVal),][1,]$common=TRUE 
  }else{    
    total[which(total$triplet==trip & total$count==maxVal),]$common=TRUE
  }    
}



total_min=total[total$common==F,]
total_max=total[total$common==T,]

chiP=c()
for (trip in unique(as.character(total_min$triplet)))
{
    p_v=chisq.test(subset(total_min,triplet==trip)$count)$p.value
    chiP=c(chiP,p_v)
}    
chiP=p.adjust(chiP,method="fdr")
total_min$Q=rep(chiP,each=2)
total_max$Q=0
total=rbind(total_min,total_max)
total=total[complete.cases(total), ] 
total$BICdiff = total$BIC2-total$BIC1
total$Qsig = total$Q<0.05
total$sig=total$BICdiff < -10
total$Order="Odonata"

total$Suborder=ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Zygoptera),1,all),"Zygoptera",
                      ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Anisoptera),1,all),"Anisoptera",
                      ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Anisozygoptera),1,any),"Anisozygoptera","RANDOM")))

total$Superfamily=replace(as.character(total$superfamily), as.character(total$superfamily)=="Anisozygoptera", "RANDOM")
total[total$Superfamily=="Cardulegastroidea","Superfamily"]="Cordulegastroidea"

total$focalclade=ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Aeshnidae),1,all),"Aeshnidae",
                      ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Gomphidae_Petaluridae),1,all),"Gomphidae+Petaluridae",
                      ifelse(apply(apply(total[,c("P1","P2","P3")],2,"%in%",Libellulidae),1,all),"Libellulidae","RANDOM")))


total$type=ifelse(total$common==TRUE & total$sig==TRUE,"Concordant",ifelse(total$common==TRUE & total$sig==FALSE,"Extreme ILS",ifelse(total$sig==FALSE & total$common==FALSE,"ILS","Introgression+ILS")))
total$Qtype=ifelse(total$common==TRUE & total$sig==TRUE,"Concordant",ifelse(total$common==TRUE & total$sig==FALSE,"Extreme ILS",ifelse(total$Qsig==FALSE & total$common==FALSE ,"ILS","Introgression+ILS")))



total_mixprop=melt(total[total$type=="Introgression+ILS" ,c("mixprop2","Order","Suborder","Superfamily","focalclade")],value.name="taxon",id=c("mixprop2"))
total_mixprop=total_mixprop[total_mixprop$taxon!="RANDOM",]
total_mixprop$variable=replace(as.character(total_mixprop$variable),as.character(total_mixprop$variable)=="focalclade","Focal clade")
total_mixprop$variable_f=factor(total_mixprop$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))


g1=ggplot(total_mixprop, aes(x=taxon, y=mixprop2))+geom_violin(fill='salmon')+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+theme(axis.text.x = element_text(size = 8,angle=10))+ylab(expression(pi[2]))+xlab("")+ggtitle("A")

total_mixprop=melt(total[total$Qtype=="Introgression+ILS" ,c("mixprop2","Order","Suborder","Superfamily","focalclade")],value.name="taxon",id=c("mixprop2"))
total_mixprop=total_mixprop[total_mixprop$taxon!="RANDOM",]
emptydata=data.frame(mixprop2=c(0,0,0,0),variable=c("Superfamily","Superfamily","focalclade","focalclade"),taxon=c("Cordulegastroidea","Lestoidea","Aeshnidae","Gomphidae+Petaluridae"))
total_mixprop=rbind(total_mixprop,emptydata)
total_mixprop$variable=replace(as.character(total_mixprop$variable),as.character(total_mixprop$variable)=="focalclade","Focal clade")
total_mixprop$variable_f=factor(total_mixprop$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))

gq1=ggplot(total_mixprop, aes(x=taxon, y=mixprop2))+geom_violin(fill='salmon')+facet_grid(~variable_f,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+theme(axis.text.x = element_text(size = 8,angle=10))+ylab(expression(pi[2]))+xlab("")+ggtitle("C")


total_p=melt(total[ ,c("type","Order","Suborder","Superfamily","focalclade")],id="type",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="focalclade","Focal clade")
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))

g2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=type))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5))+theme(axis.text.x = element_text(size = 8,angle=10),legend.position=c("top") ,legend.direction="horizontal")+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("gray48","orange","gold", "salmon"),name="")+ggtitle("B")

total_p=melt(total[ ,c("Qtype","Order","Suborder","Superfamily","focalclade")],id="Qtype",value.name="taxon")
total_p=total_p[total_p$taxon!="RANDOM",]
total_p$variable=replace(as.character(total_p$variable),as.character(total_p$variable)=="focalclade","Focal clade")
total_p$variable_f=factor(total_p$variable, levels=c("Order","Suborder","Superfamily","Focal clade"))

gq2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=Qtype))+geom_bar(position="fill")+facet_grid(~variable_f,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5))+theme(axis.text.x = element_text(size = 8,angle=10),legend.position=c("top") ,legend.direction="horizontal")+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("gray48","orange","gold", "salmon"),name="")+ggtitle("D")





quartz(width=15.4,height=10.6)
grid.arrange(g1,gq1,g2,gq2,nrow=2,ncol=2)
quartz.save("quibl_distribution.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("quibl_distribution.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)

mutate(the_rank  = rank(-value, ties.method = "random"))


#STATS for Focal clades 





Epio=total[apply(total[,c("P1","P2","P3")]=="Epiophlebia_superstes",1,any) & total$outgroup!="Epiophlebia_superstes",]
Epio$ipair=get_intropair(Epio,"Epiophlebia_superstes")

Epio_sig=Epio[Epio$sig==TRUE & Epio$common==FALSE,]


    



