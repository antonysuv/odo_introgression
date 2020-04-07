library("ape")


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



all_triplets=function(taxa_list,gene_trees)
{
    dir.create("quibl_dir")
    setwd("quibl_dir")
    taxa_combn=combn(taxa_list,m=3)
    for (i in 1:ncol(taxa_combn))
    {
        triplet=taxa_combn[,i]
        quible_trees(triplet,gene_trees,i)
    }    
    setwd("..")
}    





Lestoidea=c("Perissolestes_remotus","Synlestes_weyersii","Episynlestes_cristatus","Indolestes_peregrinus","Archilestes_grandis","Protosticta_beaumonti")
Epiophlebiidae="Epiophlebia_superstes"
Aeshnidae=c("Telephlebia_godeffroyi", "Austroaeschna_subapicalis","Gynacantha_tibiata","Anax_parthenope","Anax_walsinghami","Anax_junius","Aeshna_palmata")
Epio=c(Lestoidea,Epiophlebiidae,Aeshnidae)



all_triplets(Epio,"/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all_wboot")

quibl=read.csv("/Users/Anton/Downloads/quibl_out_all_epio_anax")

quibl$BICdiff = quibl$BIC2-quibl$BIC1