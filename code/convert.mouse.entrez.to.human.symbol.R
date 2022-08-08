
convert.mouse.entrez.to.human.symbol <- function(x){
  genesV2=read.csv("data/ATACseq/data/Human.Mouse.Rat.Ortholog.Gene.List.csv",
                   header=T)
  genesV2=genesV2[,c("Mouse_Entrez","Human_Symbol")]
  
  genesV2[,2]=toupper(genesV2[, 2])
  gene.and.CA=genesV2[match(as.numeric(x),genesV2[,1]),]
  
  print(mean(gene.and.CA[,1]==as.numeric(x),na.rm=T))
  return(gene.and.CA[,2])
}