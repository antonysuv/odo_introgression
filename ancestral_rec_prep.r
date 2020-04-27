library("ape")
library("seqinr")



#######################################################RUN IQTREE ancestral reconstruction ####################

seq=read.table("deep.state",stringsAsFactors = F)
names(seq)=c("Node","Site","State","p_A","p_C","p_G","p_T")
seq$Node=ifelse(seq$Node=="Node20",">RZ",ifelse(seq$Node=="Node66",">Outgroup",ifelse(seq$Node=="Node17",">Anisoptera",">Lestoidea")))
seq_to_fa=tapply(seq[,c("State")],seq[,c("Node")],paste,collapse="")
write.table(seq_to_fa,"deepTest.fasta",quote=F,col.names=F)