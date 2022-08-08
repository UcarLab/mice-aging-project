#setwd("~/Documents/human_aging_footprint_copy_modifiable/")

library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(ggplot2)

#source("/data/youna/human_aging_footprint/merge_TF_family.r")
source("code/merge_TF_family_human_ldk.R")
#load("/data/youna/human_aging_footprint/unexpressedgene_second_cohort.Rdata")
load("data/human_aging_footprint/unexpressedgene_second_cohort.Rdata")
cor.cutoff=0.9
Isin <- function(x,y){
    sel=!(toupper(convert.TF.name(x))%in% toupper(y))
    return(sel)
}
TFname<-function(TF){
  ### removing .RC in the TF-motif name
  for(i in 1:length(TF)){
    tmp=strsplit(TF[i],"")[[1]]
    TF[i]=paste(tmp[-(1:7)],collapse="")
  }
  return(TF)
}

combine <- function(x){
    mat=matrix(NA,nr=nrow(x),nc=4)
    tmp=c("F O","F Y","M O","M Y")
    colnames(mat)=tmp
    rownames(mat)=rownames(x)
    for(i in 1:4)
      mat[,i]=rowMeans(x[,colnames(x)==tmp[i]])
    return(mat)
  }

riskmat <- function(x,k){
  allnames=NULL

  for(i in 1:length(x)){
    load(x[i])
    allnames=c(allnames,names(result$p.bifet))
  }
  allnames=unique(allnames)
  risk=matrix(NA,nr=length(allnames)+1,nc=length(x))
#    if(k=="total.footprints") rn1="background CA" else rn1="target CA"
  rownames(risk)=c("CA",allnames)
  allsex=NULL
  for(i in 1:length(x)){
      load(x[i]) 
      allsex=c(allsex,sex)
      if(increasing) risk[1,i]=mean(reads.positive.peak) else risk[1,i]=mean(reads.negative.peak)
      
      
      if(k=="Target.proportion") risk[names(result$p.bifet),i]=result$risk[,1]
        
      if(k=="Background.proportion") risk[names(result$p.bifet),i]=result$risk[,2]
      
      if(k=="Target/Background") risk[names(result$p.bifet),i]=result$risk[,1]/result$risk[,2]
  
      if(k=="total.footprints"){
        risk[1,i]=mean(reads.background.peak)
        risk[names(result$p.bifet),i]=result$risk[,3]
      }
      if(k=="pvalue"){
        # risk[names(result$p.bifet),i] <- -log((result$p.bifet)+10^(-20))
        # added FDR
        risk[names(result$p.bifet),i] <- -log(p.adjust(result$p.bifet, "fdr")+10^(-20))
        #risk[names(result$p.bifet),i] <- -log10(p.adjust(result$p.bifet, "fdr")+10^(-20))
        
      }
  }
  colnames(risk)=allsex
  return(risk)
}


pdf(file="data/human_aging_footprint/topregulator_ldk.pdf",width=20,height=15)


for(increasing in c(TRUE,FALSE)){
    count=0
    if(increasing){
      pcutoff=0.05
      peaktype="opening peak"
    }else{
      pcutoff =0.05
      peaktype="closing peak"
    }
    matall=matall2=CA=CA2=vector("list",4)
    for(option in c("pvalue","Target.proportion","total.footprints","Target/Background")){
        count=count+1
### original young
        backgroundcutoff=0.2

        if(increasing) agegroup="O" else agegroup="Y"
        #x=dir(path="/data/youna/human_aging_footprint",pattern=paste("BiFETresultH",agegroup,sep=""),full.names=T)
        x=dir(path = "data/human_aging_footprint/originalresult",
              pattern=paste("BiFETresultH",agegroup,sep=""),full.names=T)
        #x2=dir(path="/data/youna/human_aging_footprint",pattern=paste(backgroundcutoff,".Rdata",sep=""),full.names=T)
        x2=dir(path = "data/human_aging_footprint/originalresult",
               pattern=paste(backgroundcutoff,".Rdata",sep=""),full.names=T)

        x=intersect(x,x2)
        print(x)
        young=riskmat(x,option)

        colnames(young)=paste(colnames(young),agegroup)

### opposite old

        if(increasing) agegroup="Y" else agegroup="O"
        #x=dir(path="/data/youna/human_aging_footprint",pattern=paste("BiFETresultH",agegroup,sep=""),full.names=T)
        x=dir(path = "data/human_aging_footprint/originalresult",
              pattern=paste("BiFETresultH",agegroup,sep=""),full.names=T)
        #x2=dir(path="/data/youna/human_aging_footprint",pattern="opposite.Rdata",full.names=T)
        x2=dir(path = "data/human_aging_footprint/originalresult",
               pattern="opposite.Rdata",full.names=T)

        x=intersect(x,x2)
        print(x)
        old=riskmat(x,option)
        colnames(old)=paste(colnames(old),agegroup)

                                        #young[is.na(young)]=0
                                        #old[is.na(old)]=0

        common=intersect(rownames(young),rownames(old))
        print("common")
        print(length(common))
        common=common[Isin(common,unexpressedgene)]
        print(length(common))

        
        mat=cbind(young[common,],old[common,])
        mat[is.na(mat)]=0 # why set to 0 for all? why not 0 just for p values and remove the rows which have NAs? (like in mice)
        if(option=="total.footprints"){
          totalfootprints=colSums(mat[-1,])
          save(totalfootprints,file="data/human_aging_footprint/total_footprints_no_ldk.Rdata")
        }
        # IS THIS MANUAL BONFERRONI ADJUSTMENT????
        # Assuming yes. Line will be commented and log(pcutoff) used instead
        # Above, in riskmat, p-values will be FDR adjusted
        if(option=="pvalue"){
          #mat = mat[which(rowSums(mat> -log(pcutoff/nrow(mat)))>1),]
          #mat <- mat[which(rowSums(mat> -log10(pcutoff))>1 ),]
          mat <- mat[which(rowSums(mat> -log(pcutoff))>1 ),]
          sel=rownames(mat)
        }


        matall[[count]]=mat[sel,order(colnames(mat))]
        matall2[[count]]=combine(matall[[count]])


        CA[[count]]=matall[[count]][1,]
        CA2[[count]]=matall2[[count]][1,]
        matall[[count]]=matall[[count]][-1,]
        matall2[[count]]=matall2[[count]][-1,]

        
        ## load("topFM.Rdata")
        ## topFmat=cbind(young[topF,],old[topF,])
        ## topFmat=topFmat[,order(colnames(topFmat))]
        ## topMmat=cbind(young[topM,],old[topM,])
        ## topMmat=topMmat[,order(colnames(topMmat))]
        ## commonmat=cbind(young[common,],old[common,])
        ## commonmat=commonmat[,order(colnames(commonmat))]

                                        #mat=commonmat
                                        #df=data.frame(x=1:ncol(mat),y=mat["NRF1",],type=colnames(mat))
                                        #p<-ggplot(df, aes(x=x, y=y, fill=type)) + geom_bar(stat="identity")+theme_minimal()
                                        #MAX=max(abs(c(topFmat,topMmat,commonmat)),na.rm=T);breaksList = seq(0, MAX, by = 10^(-5))


                                        #pheatmap(mat[rowSums(is.na(mat))==0,],cluster_cols = FALSE,scale="none",main=option,fontsize=5)#,color = colorRampPalette(c("blue","white","red"))(length(breaksList)),breaks = breaksList)
      }
    newgroup= merge.TF.family(matall[[2]],TFfamily)
    
    raw_names_FC <- matall2[[4]]
    raw_names_pvals <- matall2[[1]]
    save(raw_names_FC, raw_names_pvals, file = paste0("data/human_aging_footprint/raw_FC_p", increasing, ".Rdata"))
    
    for(i in 1:length(matall)){
        matall[[i]]=merge.matrix(matall[[i]],newgroup)
        matall2[[i]]=merge.matrix(matall2[[i]],newgroup)
    }
                                        #sel=c("MA00762ELK4", "MA04921JUND")    


    
    ha <- function(x,y,z){
      ano <- HeatmapAnnotation(type = colnames(x),
                               CA=z,
                               col=list(type=c("F Y"="red","F O"="blue","M Y"="green","M O"="yellow"),
                                        CA = colorRamp2(c(min(z), max(z)), c("white", "red"))),
                               show_legend=y)
      return(ano)
    }
      


    draw(Heatmap(matall[[1]],
                 cluster_columns=F,
                 heatmap_legend_param = list(title = "log BiFET pvalue"),
                 column_title="log BiFET pvalue", 
                 row_names_gp = gpar(fontsize = 9),
                 row_title=peaktype,
                 top_annotation = ha(matall[[1]],TRUE,CA[[1]])) +   
           Heatmap(matall[[2]],
                   cluster_columns=F,
                   heatmap_legend_param = list(title = "proportion of target 
                                               peaks with footprints"),
                   column_title="proportion of target peaks with footprints", 
                   row_names_gp = gpar(fontsize = 9),
                   top_annotation=ha(matall[[2]],FALSE,CA[[2]]))        +
           Heatmap(matall[[3]],
                   cluster_columns=F,
                   heatmap_legend_param = list(title = "total footprints"),
                   column_title="total footprints", 
                   row_names_gp = gpar(fontsize = 9),
                   top_annotation=ha(matall[[3]],FALSE,CA[[3]])))#+Heatmap(matall[[4]],cluster_columns=F,heatmap_legend_param = list(title = "total footprints"),column_title="total footprints", row_names_gp = gpar(fontsize = 9)))

    draw(Heatmap(matall[[1]],
                 cluster_columns=F,
                 heatmap_legend_param = list(title = "log BiFET pvalue"),
                 column_title="log BiFET pvalue", 
                 row_names_gp = gpar(fontsize = 9),
                 row_title=peaktype,
                 top_annotation = ha(matall[[1]],TRUE,CA[[1]]),
                 show_row_names=TRUE)  +
           Heatmap(matall2[[1]],
                   cluster_columns=F,
                   heatmap_legend_param = list(title = "combined log BiFET pvalue"),
                   column_title="combined log BiFET pvalue", 
                   row_names_gp = gpar(fontsize = 9),
                   row_title=peaktype,
                   top_annotation = ha(matall2[[1]],FALSE,CA2[[1]]))    )


    draw(Heatmap(matall[[4]],
                 cluster_columns=F,
                 heatmap_legend_param = list(title = "Target/Background FC"),
                 column_title="Target/Background FC", 
                 row_names_gp = gpar(fontsize = 9),
                 row_title=peaktype,
                 top_annotation = ha(matall[[4]],TRUE,CA[[4]]),
                 show_row_names=TRUE)  +
           Heatmap(matall2[[4]],
                   cluster_columns=F,
                   heatmap_legend_param = list(title = "combined FC"),
                   column_title="combined FC", 
                   row_names_gp = gpar(fontsize = 9),
                   row_title=peaktype,
                   top_annotation = ha(matall2[[4]],FALSE,CA2[[4]]))    )
    
    combined_p_vals <- matall2[[1]]
    combined_FC <- matall2[[4]]
    CA_p <- CA2[[1]]
    CA_FC <- CA2[[4]]
    
    save(combined_p_vals, combined_FC, CA_p, CA_FC, ha,  
         file = paste0("data/human_aging_footprint/combined_human_FC_p_CA_ha", increasing, ".Rdata"))
    
}
dev.off()


# load("data/human_aging_footprint/samplereadsinfo.Rdata")
# mean(as.numeric(sample[,3])[sample[,2]=="F young"])
# mean(as.numeric(sample[,3])[sample[,2]=="M young"])
# mean(as.numeric(sample[,3])[sample[,2]=="F old"])
# mean(as.numeric(sample[,3])[sample[,2]=="M old"])





