library("ape")
library('ggplot2')
library('pals')
library('reshape2')
library('dplyr')
library('gridExtra')
library("svMisc")


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
        progress(i,ncol(taxa_combn))
        triplet=taxa_combn[,i]
        quible_trees(triplet,gene_trees,i)
    }    
    setwd("..")
}  

all_triplets_az=function(taxa_list,gene_trees,dir_name)
{
    dir.create(dir_name)
    setwd(dir_name)
    taxa_combn=combn(taxa_list,m=2)
    taxa_combn=rbind(taxa_combn,rep("Epiophlebia_superstes",ncol(taxa_combn)))
    for (i in 1:ncol(taxa_combn))
    {
        progress(i,ncol(taxa_combn))
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

Anisozygoptera=drop.tip(tt,c("Epiophlebia_superstes","Ephemera_danica", "Isonychia_kiangsinensis"))$tip.label
Zygoptera=c(extract.clade(tt,88)$tip.label)
Anisoptera=c(extract.clade(tt,136)$tip.label)

all_triplets_az(Anisozygoptera,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Anisozygoptera_quibl")
all_triplets(Zygoptera,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Zygoptera_quibl")
all_triplets(Anisoptera,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot","Anisoptera_quibl")



################################################################ RUN QuiBL ########################################################################