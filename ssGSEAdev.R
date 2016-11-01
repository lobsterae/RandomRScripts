library(GSVA)
library(pheatmap)
library(RankProd)
library(ggplot2)
library(ggbio)
library(BiocParallel)
library(plotly)
#------------------------------------------------------
n.procs=12
register(MulticoreParam(workers=n.procs))
#------------------     READ ARGS  --------------------
options(echo=TRUE)
#------------------------------------------------------
args=commandArgs(trailingOnly=TRUE)
#------------------------------------------------------
print(args)
#------------------------------------------------------
gmt.file =args[1]
expr.file=args[2]
ctrl.file=args[3]
#------------------------ LOAD DATA ---------------------
#--------------------------------------------------------
getEnrichment <- function(eset, gmt, 
			  weighted.score.type = 1, 
			  rnaseq=True, method="ssgsea", 
			  correl.vector = NULL)
		 {
			rownames(eset)=tolower(rownames(eset))
			scores=NULL
			scores=apply(eset,2,function(g){
							gsva(as.vector(g),gmt,
					     		rnaseq,
							    gmt,
				      		method=method)
						      })
			return(scores)
     }  
#-----------------------------------------------------------
scoreEnrichment <- function(scores.test,scores.ctrl)
		   {
			stats=list()
				for( i in 1:nrow(scores.test))
				{
				  for(j in 1:ncol(scores.test))
				  {
				   stats[[i]][[j]]$WCOX=wilcox.test(scores.test[i,j],scores.ctrl[i,],
								exact=TRUE,alt="greater")
				   stats[[i]][[j]]$Z$score=abs(scores.test[i,j]-mean(scores.ctrl[i,]))/sd(scores.ctrl[i]))
					 stats[[i]][[j]]$Z$pval=1-pnorm(abs(stats[[i]][[j]]$Z$score))
				  }
				}

      	return(lapply(stats,
      	              lapply(stats,
      	                     function(s){
      	                       return(c(s$WCOX$statistic,
      	                                s$WCOX$p.value,
      	                                s$Z$score,
      	                                s$Z$pval))}
      	                     )))

		   }
#------------------------------------------------------------
loadGMT<-function(gmtfile)
	 {
		gmt=NULL
		gmt=readLines(gmt)
		gmt=lapply(gmt,function(x){return(unlist(strsplit(x)))})
		nms=lapply(gmt,function(x){return(x[1])})
		gmt=lapply(gmt,function(x){return(tolower(x[2:length(x)]))})

		names(gmt)=nms
		return(gmt)	
	 }
#-------------------------------------------------------------
loadExpr<-function(exprfile)
 	  {
		eset=NULL
		eset=read.table(exprfile,header=T,row.names=1,sep="\t")
		return(eset)
	  }

#-------------------------------------------------------------
getGeneRanks<-function(data=expr,cl=rep(1,ncol(expr)),logged=T,na.rm=F,huge=F)
{
		exprRP=apply(expr,2,function(x){
			return(RP(data=x,cl=cl,logged=logged,na.rm=na.rm,huge=huge))
		        })

		exprRP=lapply(exprRP,as.data.frame)	
		return(exprRP)
}
#---------------------------------------------------------
printResults<-function(expr,scores,gmt,pthresh)
{
	#scores.signif=return(scores$)
}
#---------------------------------------------------------
writeOutput<-function(test,ctrl,gsva.test,gsva.ctrl,
			scores.test,scores.ctrl,gmt,pthresh)
	   {
		# scores that are signif 
		# signif.samples=lapply(scores,function(x){return(which(x[,2]<0.05))})
		
		NAMES.GMT=names(gmt)
		SIGNIF.SMPLS=lapply(scores.test,
				function(x){lapply(x,
					    function(y){return(length(which(y[,2]<pthresh)))})}) # number of gmt hits per sample
	
		SIGNIF.SMPLS.IDX=which(unlist(SIGNIF.SMPLS)>0)		
		print("No. of samples with hits : "+length(SIGNIF.SMPLS.IDX))
		scores.test=scores.test[SIGNIF.SMPLS.IDX]
		gsva.test=gsva.test[,SIGNIF.SMPLS.IDX]
		test=test[,SIGNIF.SMPLS.IDX]
		NAMES.SMPLS=colnames(test)

		for( i in 1:length(scores.test))
		{
		  #one list for each smpl
		  smpl.name=NAMES.SMPLS[i]
		  
		  signif.gmt=which(scores[[i]][,2]< pthresh)
		  names.gmt=NAMES.GMT[signif.gmt]
		  nhits=length(names.gmt)
		  
		  if(length(signif.gmt)>0)
      {
			 for( j in 1:length(signif.gmt))
			 {
			 
			  write.table( c(smpl.name,
				 	names.gmt[j],
				 	length(gmt[signif.gmt[j]]),
				 	gsva.test[i,j],
				 	scores.test[[i]][signif.gmt[j],2],
				 	scores.test[[i]][signif.gmt[j],1]
					                                  ),
					file="SAMPLE_HITS.txt",
					append=TRUE,
					header=FALSE)
        }
      } 	
  }	
		
		
}
#------------------------------------------------------
#gmt      =loadGMT(gmt.file)
#test.expr=loadExpr(expr.file)
#ctrl.expr=loadExpr(ctrl.file)
#------------------------------------------------------
#
#----------------      MAKE PLOTS         -------------
#GSVA.SCORES=getEnrichment(test.expr,gmt)
#GSVA.CTRL  =getEnrichment(ctrl.expr,gmt)
#STATS.TEST =scoreEnrichment(GSVA.SCORES,GSVA.CTRL)
#GENE.RANKS.EXPR =getGeneRanks(test.expr)
#GENE.RANKS.CTRL =getGeneRanks(test.ctrl)
#GENE.PLOTS      =writeOutput( test.expr,test.ctrl,
#                                              GENE.RANKS.EXPR,
#                                                GENE.RANKS.CTRL,
#                                                GSVA.SCORES,GSVA.CTRL,
#                                                gmt,pthresh)
#--------------------------------------------------------
