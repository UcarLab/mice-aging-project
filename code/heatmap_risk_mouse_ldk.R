library(ComplexHeatmap)
library(circlize)
source("code/merge_TF_family_human_ldk.R")
load("data/mouse_aging_footprint/mouse_unexpressed_genelist.Rdata")
cor.cutoff=0.9

MERGE=FALSE
#MERGE <- TRUE
get.strain <- function(x){
  y=x
  for(i in 1:length(x)){
      tmp <- strsplit(x[i],"_merged")[[1]][1]
      tmp2=strsplit(tmp,"_")[[1]]
      y[i]=tmp2[length(tmp2)]
    }
  return(y)
}

remove.space <- function(x){
  for(i in 1:length(x)){
    x[i]=paste(strsplit(x[i],"var2")[[1]][1])

    tmp=strsplit(x[i]," ")[[1]]
    tmp[tmp=="FOSJUN"]="JUNFOS"
    tmp[tmp=="Tal1Gata1"]="GATA1TAL1"
    
    x[i]=paste(tmp,collapse="")
  }
  return(x)
}


Isin <- function(x,y){
  sel=!(toupper(convert.TF.name(x))%in% toupper(y))
  return(sel)
}

TFname<-function(TF) ### removing .RC in the TF-motif name
{
  if(length(TF)>0)
    for(i in 1:length(TF))
      {
        tmp=strsplit(TF[i],"")[[1]]
        TF[i]=paste(tmp[-(1:7)],collapse="")
      }
  return(TF)
}


without.sex <- function(x){
  y=x
  for(i in 1:length(x))
    y[i]=paste(strsplit(x[i]," ")[[1]][-1],collapse=" ")
      
  return(y)
}

combine <- function(x){
  mat=matrix(NA,nr=nrow(x),nc=4)
  colnames(mat)=c("B6 3mo","B6 18mo","NZO 3mo","NZO 18mo")
  
  rownames(mat)=rownames(x)

  mat[,"NZO 3mo"]=(x[, "F NZO 3mo"]+x[, "M NZO 3mo"])/2
  mat[,"B6 3mo"]=(x[, "F B6 3mo"]+x[, "M B6 3mo"])/2
  

  mat[,"NZO 18mo"]=(x[, "F NZO 18mo"]+x[, "M NZO 18mo"])/2
  mat[,"B6 18mo"]=(x[, "F B6 18mo"]+x[, "M B6 18mo"])/2
  
  ## mat[,"F/M NZO"]=(x[, "F NZO 18mo"]/x[, "F NZO 3mo"]+x[, "M NZO 18mo"]/x[, "M NZO 3mo"])/2
  ## mat[,"F/M B6"]=(x[, "F B6 18mo"]/x[, "F B6 3mo"]+x[, "M B6 18mo"]/x[, "M B6 3mo"])/2
  
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
    
    allsex=c(allsex,paste(sex,get.strain(x[i]),collapse=" "))
    
    if(increasing) risk[1,i]=mean(reads.positive.peak) else risk[1,i]=mean(reads.negative.peak)
    
    
    if(k=="Target.proportion")  risk[names(result$p.bifet),i]=result$risk[,1]
      
    if(k=="Background.proportion")    risk[names(result$p.bifet),i]=result$risk[,2]
    
    if(k=="Target/Background")    risk[names(result$p.bifet),i]=result$risk[,1]/result$risk[,2]

    if(k=="total.footprints"){
      risk[1,i]=mean(reads.background.peak)
      risk[names(result$p.bifet),i]=result$risk[,3]
    }
    if(k=="pvalue")    risk[names(result$p.bifet),i] <- 
      -log(p.adjust(result$p.bifet,"fdr")+10^(-20)) # -log10(p.adjust(result$p.bifet,"fdr")+10^(-20))
  }
  colnames(risk)=allsex
  rownames(risk)=remove.space(rownames(risk))
  return(risk)
}

library(pheatmap)

library(ggplot2)

for(tissue in c("spleen","PBL","memory","naive")){

  pdf(file=paste("data/mouse_aging_footprint/topregulator",tissue,"_ldk.pdf",sep=""),width=20,height=15)

  for(increasing in c(TRUE,FALSE)){
    count=0
    if(increasing){
      pcutoff=0.05
      peaktype="opening_peak"
    }else{
      pcutoff =0.05
      peaktype="closing_peak"
    }
    
    matall=matall2=CA=CA2=vector("list",4)
    
    for(option in c("pvalue","Target.proportion","total.footprints","Target/Background")){
        count=count+1
### original young
        backgroundcutoff=0.2

        if(increasing) agegroup="18mo" else agegroup="3mo"

        x=dir(path = "data/mouse_aging_footprint/F5_inputs/",
              pattern=paste("BiFETresultsub_",tissue,"_",agegroup,sep=""),
              full.names=T)

        x2=dir(path = "data/mouse_aging_footprint/F5_inputs/", 
               pattern=paste(backgroundcutoff,".Rdata",sep=""),
               full.names=T)

        x=intersect(x,x2)
        print(x)
        young=riskmat(x,option)

        colnames(young)=paste(colnames(young),agegroup)

### opposite old

        if(increasing) agegroup="3mo" else agegroup="18mo"
        x=dir(path = "data/mouse_aging_footprint/F5_inputs/", 
              pattern=paste("BiFETresultsub_",tissue,"_",agegroup,sep=""),
              full.names=T)
        
        x2=dir(path = "data/mouse_aging_footprint/F5_inputs/",
               pattern="opposite.Rdata",
               full.names=T)

        x=intersect(x,x2)
        print(x)
        old=riskmat(x,option)
        colnames(old)=paste(colnames(old),agegroup)

        #young[is.na(young)]=0
        #old[is.na(old)]=0

        common=intersect(rownames(young),rownames(old))
        print("common")
        print(length(common))
        common=common[Isin(common,unexpgene[[tissue]])]
        print(length(common))
        
        mat=cbind(young[common,],old[common,])
        
        # if p value - changes NAs to zeros
        # if any other - removes rows which have an NA
        
        ###This should be changed to get rid of###
        if(option=="pvalue") mat[is.na(mat)]=0 else mat=mat[rowSums(is.na(mat))==0,]
        
        if(option=="Target/Background"){
            FC.pvalue=rep(1,nrow(mat))
            for(fci in 1:nrow(mat)){
              tmp=try(t.test(log(mat[fci,1:ncol(young)]),log(mat[fci,-(1:ncol(young))]),paired=T)$p.value,TRUE)
              #tmp=try(t.test(log10(mat[fci,1:ncol(young)]),log10(mat[fci,-(1:ncol(young))]),paired=T)$p.value,TRUE)
              if(!inherits(tmp, "try-error")) FC.pvalue[fci]=tmp
            }
            selFC=rownames(mat)[FC.pvalue<1.1]
        }

        if(option=="total.footprints"){
          totalfootprints=colSums(mat[-1,])
          save(totalfootprints,
               file="data/mouse_aging_footprint/total_footprints_no_ldk.Rdata")
        }
          

        if(option=="pvalue"){
          mat=rbind(mat[which(rowSums(mat> -log(pcutoff))>1 ),])
          #mat=rbind(mat[which(rowSums(mat> -log10(pcutoff))>1 ),])
          selBiFET=rownames(mat)
          if(nrow(mat)>1){
            write.csv(rbind(exp(-mat[-1,])),
                      file=paste("data/mouse_aging_footprint/",
                                 tissue,peaktype,"regulator_BiFET_FDR_ldk.csv",
                                 sep=""),
                      quote=F)  
          } 
        }
        
        matall[[count]] = mat[,c("F B6 3mo","M B6 3mo","F B6 18mo","M B6 18mo", "F NZO 3mo",  "M NZO 3mo","F NZO 18mo","M NZO 18mo"  )]
        
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
      
      names(matall) <- c("pvalue","Target.proportion","total.footprints","Target/Background")
      
      sel=intersect(selFC,selBiFET)
      save(sel,matall,
           file=paste("data/mouse_aging_footprint/matall",increasing,"_ldk.Rdata",sep="_"))
      write.csv(toupper(TFname(sel[sel!="CA"])),
                file=paste("data/mouse_aging_footprint/",
                           tissue,peaktype,"regulator_ldk.csv",sep=""),
                quote=F,row.names=F)
    }


    matall2=CA=CA2=vector("list",4)
    combinedheatmap=vector("list",2)
    names(combinedheatmap)=c("increasing","decreasing")
    for(increasing in c(TRUE,FALSE)){
      
      if(increasing){
        peaktype="opening_peak"
      }else{
        peaktype="closing_peak"
      }
  
      load(paste("data/mouse_aging_footprint/matall",increasing,"_ldk.Rdata",sep="_"))
      if(length(sel)>0){
        sel=unique(c("CA",sel))
  
        for(count in 1:length(matall)){
          matall[[count]]=matall[[count]][sel,]
          matall2[[count]]=combine(matall[[count]])

          tmprowname=rownames(matall[[count]])
      
          CA[[count]]=matall[[count]][1,]
          CA2[[count]]=matall2[[count]][1,]
          matall[[count]]=rbind(matall[[count]][-1,])
          matall2[[count]]=rbind(matall2[[count]][-1,])
          rownames(matall[[count]])=tmprowname[-1]
          rownames(matall2[[count]])=tmprowname[-1]
      
        }
        
        names(matall2) <- paste0("combined_", names(matall))
        
        if(MERGE){
          newgroup= merge.TF.family(matall[[2]],TFfamily)
    
          if(!is.null(newgroup))      
            for(i in 1:length(matall)){
              matall[[i]]=merge.matrix(matall[[i]],newgroup)
              matall2[[i]]=merge.matrix(matall2[[i]],newgroup)
            }
          }
         #sel=c("MA00762ELK4", "MA04921JUND")    
    
      ha <- function(x,y,z){
        heat_annot <- HeatmapAnnotation(type = without.sex(colnames(x)),
                                    CA=z,
                                    col=list(type=c("B6 3mo"="red",
                                                    "B6 18mo"="blue",
                                                    "NZO 3mo"="green",
                                                    "NZO 18mo"="yellow"),
                                             CA= colorRamp2(c(min(z), 
                                                              max(z)), 
                                                            c("white", "red"))),
                                    show_legend=y)
        return(heat_annot)
      }
         
           
      draw(Heatmap(matall[[1]],cluster_columns=F,
             heatmap_legend_param = list(title = "log BiFET pvalue"),
             column_title="log BiFET pvalue", 
             row_names_gp = gpar(fontsize = 9),
             row_title=peaktype, 
             show_row_names=TRUE,
             top_annotation = ha(matall[[1]],TRUE,CA[[1]])) +   
           Heatmap(matall[[2]],cluster_columns=F,
             heatmap_legend_param = 
             list(title = "proportion of target peaks with footprints"),
             column_title="proportion of target peaks with footprints", 
             row_names_gp = gpar(fontsize = 9),
             show_row_names=TRUE,
             top_annotation=ha(matall[[2]],FALSE,CA[[2]])) +
           Heatmap(matall[[3]], cluster_columns=F,
              heatmap_legend_param = list(title = "total footprints"),
              column_title="total footprints", 
              row_names_gp = gpar(fontsize = 9),
              show_row_names=TRUE,
              top_annotation=ha(matall[[3]],FALSE,CA[[3]])),
           auto_adjust = FALSE)#+Heatmap(matall[[4]],cluster_columns=F,heatmap_legend_param = list(title = "total footprints"),column_title="total footprints", row_names_gp = gpar(fontsize = 9)))

    draw(Heatmap(matall[[1]],cluster_columns=F,
                 heatmap_legend_param = list(title = "log BiFET pvalue"),
                 column_title="log BiFET pvalue", 
                 row_names_gp = gpar(fontsize = 9),
                 row_title=peaktype,
                 top_annotation = ha(matall[[1]],TRUE,CA[[1]]),
                 show_row_names=TRUE)  +
           Heatmap(matall2[[1]], cluster_columns=F, 
                   heatmap_legend_param = list(title = "combined log BiFET pvalue"),
                   column_title="combined log BiFET pvalue",
                   row_names_gp = gpar(fontsize = 9),
                   #top_annotation = ha(matall2[[1]],FALSE,CA2[[1]]),
                   row_title=peaktype),
         auto_adjust = FALSE)


      draw(Heatmap(matall[[4]], cluster_columns=F,
             heatmap_legend_param = list(title = "Target/Background FC"),
             column_title="Target/Background FC", 
             row_names_gp = gpar(fontsize = 9),
             row_title=peaktype,top_annotation = ha(matall[[4]],TRUE,CA[[4]]),
             show_row_names=TRUE)  +
           Heatmap(matall2[[4]],cluster_columns=F,
             heatmap_legend_param = list(title = "combined FC"),
             column_title="combined FC", 
             row_names_gp = gpar(fontsize = 9),row_title=peaktype),
           auto_adjust = FALSE)#,top_annotation = ha(matall2[[4]],FALSE,CA2[[4]]))    )
      
      save(matall, matall2, file = paste0("data/mouse_aging_footprint/",
                                          tissue, "_", increasing,"_", 
                                          "_matall_lists_ldk.Rdata"))
      }
      
    if(increasing){
      combinedheatmap[["increasing"]]=matall2[[4]]
    }else{
      combinedheatmap[["decreasing"]]=matall2[[4]]  
    }
  }
  dev.off()
  
  save(combinedheatmap,file=paste("data/mouse_aging_footprint/topregulator",tissue,"_ldk.Rdata",sep=""))
}
