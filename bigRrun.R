#--- big R Run -------------------------
# make the grand results and stats table
#---------------------------------------
library(stats)
library(DOSE)
#---------------------------------------
args = options(trailing=TRUE)
gmt  = args[1]
name = args[2]
#---------------------------------------
loadGMT(gmt)
load(paste(name,".Rda",sep=""))
#---------------------------------------
PVAL=0.05
SIGNIF.TOT=length(unique(RP[[1]][which(RP[[5]]<0.05)]))
STATS=list()

map2geneLen=function(e.gene)
{
  library(org.Mm.eg.db)
  tmp=as.list(org.Mm.egACCNUM)
  keys=mappedKeys(org.Mm.eg.db)
}

getNES=function(gene.list,gene.set,weighted=1,correl.vector=NULL,permut=100)
{
  
  
}


getS2N=function(gene.list,gene.set,weighted=1,correl.vector=NULL,permut=100)
{
  
  
}

#--------------------------------------------------------------------
gseaScores=function (geneList, geneSet, exponent = 1)
{
  geneSet      <- intersect(geneSet, names(geneList))
  N            <- length(geneList)
  Nh           <- length(geneSet)
  Phit         <- Pmiss <- numeric(N)
  hits         <- names(geneList) %in% geneSet
  Phit[hits]   <- abs(geneList[hits])^exponent
  NR           <- sum(Phit)
  Phit         <- cumsum(Phit/NR)
  Pmiss[!hits] <- 1/(N - Nh)
  Pmiss        <- cumsum(Pmiss)
  runningES    <- Phit - Pmiss
  max.ES       <- max(runningES)
  min.ES       <- min(runningES)
  
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  }
  else {
    ES <- min.ES
  }

  df <- data.frame(x = seq_along(runningES), runningScore = runningES,
                 positionps - = as.integer(hits))
  return(df)

}
#---- function getScores using full score distro
getScores=function(exprs,gmt,nms)
{
  tmp=exprs[,1]
  names(tmp)=nms
  tmp=tmp[order(tmp,decreasing=TRUE)]  
  m=match(gmt,nms)
  scores=NULL
  
  for(i in 1:(length(tmp)-length(m)))
  {
    idx =seq(1,length(tmp),by=1)
    swap=seq(i,(i+length(m))-1,by=1)
    idx[swap]= m
    idx[m]=swap
    
    scores=cbind(scores,
             gsva(expr=as.matrix(tmp[idx]),
             gset.idx.list=list(gmt=lapply(gmt,tolower)),method="ssgsea"))
    i=i+50
  }
  
  return(scores)
  
}

for(i in 1:length(RP))
{
  STATS[[i]]=list()
  
  for(j in 1:length(gmt))
  {
    STATS[[i]][[j]]=list()
    
    m=match(gmt.filt[[j]],dRP.genes)
    
    STATS[[i]][[j]]$RP   = dRP[[i]][na.omit(m),]
    STATS[[i]][[j]]$LEN  = length(na.omit(m))
    STATS[[i]][[j]]$RAWL = length(m)
    STATS[[i]][[j]]$PROP = length(which(STATS[[i]][[j]]$RP[,5]<0.05))
    
    a=STATS[[i]][[j]]$PROP
    b=length(na.omit(m))-a
    c=SIGNIF.TOT-a
    d=(length(dRP[[i]][,1])-SIGNIF.TOT)-c
    
    STATS[[i]][[j]]$CHISQ= fisher.test(rbind(c(a,b),c(c,d)),alt="greater")
    STATS[[i]][[j]]$PVAL = STATS[[i]][[j]]$CHISQ$p.value
    STATS[[i]][[j]]$LOD  = STATS[[i]][[j]]$CHISQ$estimate
    STATS[[i]][[j]]$BONF = STATS[[i]][[j]]$CHISQ$p.value*length(gmt)
    STATS[[i]][[j]]$BH   = STATS[[i]][[j]]$CHISQ$p.value
  }
  
}

N=1

for(i in 1:length(STATS))
{
  
  for(j in 1:length(gmt.filt))
  {
    gmtname=gmt.names[1]
    
    write.table(cbind(gmtname,STATS[[i]][[j]]$RP),
                file="GMT_RESULTS_ALL.txt",append=TRUE,row.names=FALSE)
    ans="no"
    
    
    GSEA.obj=new("gseaResult",
                 result      =data.frame(ID="test",
                                         Description    ="GSEA.data",
                                         setSize        =length(gmt[[j]]),
                                         enrichmentScore=GSEA[i,j],
                                         NES            =GSEAnorm[i,j],
                                         pvalue         =STATS[[i]][[j]]$p.value,
                                         p.adjust       =STATS[[i]][[j]]$BONF,
                                         qvalues        =STATS[[i]][[j]]$BONF),
                 setType     ="custom",
                 geneList    =STATS[[i]][[j]]$RP[which(STATS[[i]][[j]]$RP[,7]<0.05),],
              #   permScores  =getScores(STATS[[i]][[j]]$RP),
                 params      =list(minGSSize    =5,
                                   exponent     =1e+5,
                                   pAdjustMethod="fdr",
                                   nPerm        ="100",
                                   pvalueCutoff ="0.05",
                                   organism     ="mm10",
                                   setType      ="custom"),
                 geneSets=gmt)
    
    if(STATS[[i]][[j]]$PVAL<0.05)
    {
      ans="yes"
      gseaplot(GSEA.obj,N,by="all")
      N=N+1
    }
    
    write.table(c(gmtname,
                  STATS[[i]][[j]]$LOD,
                  STATS[[i]][[j]]$PVAL,
                  STATS[[i]][[j]]$BONF,
                  STATS[[i]][[j]]$PROP,
                  STATS[[i]][[j]]$LEN,
                  STATS[[i]][[j]]$RP[,9],
                  GSVA_SCORES[i][j],ans),
                  file     ="GMT_SUMMARY_RESULTS.txt",
                  append   =TRUE,
                  row.names=FALSE)
  }
  
  
}
