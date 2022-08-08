unexpgene=vector("list",4)
names(unexpgene)=c("memory","naive","PBL","spleen")

unexpgene[["memory"]]= unique(convert.mouse.entrez.to.human.symbol(annotation[,1][rowSums(bed[,TISSUE=="CD8.memory"]>1)<1]))
unexpgene[["naive"]]= unique(convert.mouse.entrez.to.human.symbol(annotation[,1][rowSums(bed[,TISSUE=="CD8.naive"]>1)<1]))
unexpgene[["PBL"]]= unique(convert.mouse.entrez.to.human.symbol(annotation[,1][rowSums(bed[,TISSUE=="PBL"]>1)<1]))
unexpgene[["spleen"]]= unique(convert.mouse.entrez.to.human.symbol(annotation[,1][rowSums(bed[,TISSUE=="spleen"]>1)<1]))

save(unexpgene,file="data/RNAseq/results/mouse_unexpressed_genelist.Rdata")
