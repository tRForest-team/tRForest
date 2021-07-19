library(dplyr)
library(BiocManager)
library(clusterProfiler)
library(speciesorgdb) 
library(stringr)
library(cowplot)
library(ggplot2)

trforestplot=function(filename) {
  
  trforestinput=read.csv(file = filename)
  names(trforestinput)[names(trforestinput)=="name"]="gene_name"
  my.symbols=(trforestinput$gene_name)
  my.symbols=my.symbols[!duplicated(my.symbols)]
  
  if(length(my.symbols)>5)
  {
    entrezoutput=bitr(my.symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=speciesorgdb)
    names(entrezoutput)[names(entrezoutput)=="SYMBOL"]="gene_name"
    names(entrezoutput)[names(entrezoutput)=="ENTREZID"]="entrez_id"
    
    fulloutput=merge(trforestinput, entrezoutput, by="gene_name")
    fulloutput=distinct(fulloutput)
    entrezvector=fulloutput[,24]
    entrezvector=na.omit(entrezvector)
    entrezvector=unique(entrezvector)
    
    gooutputbp=enrichGO(entrezvector,speciesorgdb,keyType="ENTREZID", ont="BP", qvalueCutoff = 1, pvalueCutoff = 1)
    dotplotbp=dotplot(gooutputbp, showCategory=10)
    
    gooutputmf=enrichGO(entrezvector,speciesorgdb,keyType="ENTREZID", ont="MF", qvalueCutoff = 1, pvalueCutoff = 1)
    dotplotmf=dotplot(gooutputmf, showCategory=10)
    
    gooutputcc=enrichGO(entrezvector,speciesorgdb,keyType="ENTREZID", ont="CC", qvalueCutoff = 1, pvalueCutoff = 1)
    dotplotcc=dotplot(gooutputcc, showCategory=10)
    
    gsgooutputbp=setReadable(gooutputbp, OrgDb=speciesorgdb, 'ENTREZID')
    goplotbp=cnetplot(gsgooutputbp)
    
    gsgooutputmf=setReadable(gooutputmf, OrgDb=speciesorgdb, 'ENTREZID')
    goplotmf=cnetplot(gsgooutputmf)
    
    gsgooutputcc=setReadable(gooutputcc, OrgDb=speciesorgdb, 'ENTREZID')
    goplotcc=cnetplot(gsgooutputcc)
    
    allplotbp=plot_grid(dotplotbp, goplotbp)+theme(plot.margin=margin(50,10,50,10))
    allplotmf=plot_grid(dotplotmf, goplotmf)+theme(plot.margin=margin(50,10,50,10))
    allplotcc=plot_grid(dotplotcc, goplotcc)+theme(plot.margin=margin(50,10,50,10))
    
    trfnamepattern<-("profiles_....")
    trfname=str_extract(filename,trfnamepattern)
    trfname=sub('........_', '', trfname)
    trfname=sub("^", "Species trf-", trfname)
    trfname=paste0(trfname, " Ontology Analysis")
    
    trflabel<- ggdraw() + 
      draw_label(trfname, fontface = 'bold', x = 0, hjust = 0, size=30) +
      theme(plot.margin = margin(50,10,50,10))
    
    plots=align_plots(allplotbp, allplotmf, allplotcc, align='v', axis='l')
    
    threeplots=plot_grid(allplotbp, allplotmf, allplotcc, labels=c('Biological Process', 'Molecular Function', 'Cellular Component'), label_size=15, ncol=1, nrow=3)
    
    trfplot=plot_grid(trflabel, threeplots, ncol=1, rel_heights=c(0.05, 1))
    
    trfname2=str_extract(filename,trfnamepattern)
    trfname2=sub('........_', '', trfname2)
    trfname2=sub("^", "/Users/lainemarrah/Desktop/speciestRForestPlots/trf-", trfname2)
    trfname2=paste0(trfname2, "_plot.png")
    
    
    ggsave(trfname2, plot=trfplot, width=50, height=70, units="cm")
  }
  
  else {
    trfnamepattern<-("profiles_....")
    trfname=str_extract(filename,trfnamepattern)
    trfname=sub('........_', '', trfname)
    trfname=sub("^", "Insufficient genes for trf-", trfname)
    print(trfname)
  }
  
}

trforestspecies=list.files(path="/Users/lainemarrah/Desktop/speciesfiles", full.names=T)

for (f in trforestspecies){
  trforestplot(f)
}






