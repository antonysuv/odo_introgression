# odo_phylogeny_rnaseq

MCMC.tree.plot(zz,analysis.type = "MCMCtree",plot.type = "phylogram", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "purple", alpha.f = 0.7), no.margin = T, cex.labels = 0.01,cex.tips = 0.6,ladderize.tree = F,pos.age=-7,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6)

# Commands

java -jar ~/soft/phylonet/PhyloNet_3.8.0.jar summary.nex

python ~/soft/HyDe/scripts/run_hyde_mp.py -i ../SuperMatrix_50BUSCO_dna_pasta_ali_trim.phy -m ../taxa_map.txt -o P84_Ephemera_danica -n 85 -t 85 -s 2167861 -j 24  --prefix hyde_all_tri

python2.7 ~/soft/QuIBL/QuIBL.py ./sampleInputFile.txt

