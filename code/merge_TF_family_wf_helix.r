convert.TF.name <- function(x){
  for(i in 1:length(x)){
    x[i]=paste(strsplit(x[i],"")[[1]][-(1:7)],collapse="")
  }
  return((x))
}

ISfamily <- function(x,y,z){
    x=paste(strsplit(x,"")[[1]][-(1:7)],collapse="")
    y=paste(strsplit(y,"")[[1]][-(1:7)],collapse="")
    
    res=FALSE
    for(i in 1:length(z)){
      tmp=strsplit(z[i],",")[[1]]
      if(sum(toupper(tmp)==toupper(x))>0 & sum(toupper(tmp)==toupper(y))>0){
        res=TRUE  
      }
    }
    if(toupper(x)==toupper(y)) res=TRUE
    return(res)
  }

#TFfamily=as.matrix(read.delim("/data/youna/human_aging_footprint/motif_families_Feb2019.txt",header=T))[,"TFNames"]

TFfamily=as.matrix(read.delim("data/human_aging_footprint/motif_families_Feb2019.txt",header=T))[,"TFNames"]

merge.TF.family <- function(x,y){

    COR=cor(t(x))
    COR[upper.tri(COR,diag=T)]=0
    pair=which(COR>cor.cutoff,arr.ind=T)
    if(nrow(pair)==0){
        return(NULL)
    }else{
    membership=NULL
    for(i in 1:nrow(pair)){
      if(ISfamily(rownames(x)[pair[i,1]],rownames(x)[pair[i,2]],y)) membership=rbind(membership,pair[i,])  
    }
    
    if(is.null(membership)){
      return(NULL)
    }else{
      count=0
      group=vector("list",1)
      while(nrow(membership)>0){
        count=count+1
        sel1=which(membership==membership[1,1],arr.ind=T)[,1]
        sel2=which(membership==membership[1,2],arr.ind=T)[,1]

        tmp=sort(unique(c(sel1,sel2)))

        group[[count]]=unique(as.vector(membership[tmp,]))

        membership=rbind(membership[-tmp,])
      }
    
      newgroup=vector("list",1)
      count=0

      group0=group
      while(length(group)>0){
        count=count+1
        tmp=group[[1]]
          #print(tmp)
        group[[1]]=NULL
        i=1
        while(i <=length(group)){
          if(length(intersect(tmp,group[[i]]))>0){
            tmp=c(tmp,group[[i]])
            group[[i]]=NULL
          }else{
            i=i+1
          } 
        }
        newgroup[[count]]=unique(tmp)
      }
      return(newgroup)
    }
  }
}

merge.matrix <- function(x,newgroup){
  rownames(x)=convert.TF.name(rownames(x))

  if(is.null(newgroup)){
    newx=x
  }else{
    newx=NULL
    for(i in 1:length(newgroup)){
      newx=rbind(newx,colMeans(x[newgroup[[i]],]))
      rownames(newx)[i]=paste(sort(unique(rownames(x)[newgroup[[i]]])),collapse=",")
    }
    newx=rbind(newx,x[-unique(unlist(newgroup)),])
  }
  return(newx)
}
  
