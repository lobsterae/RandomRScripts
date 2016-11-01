#-------------------------------------------
library(timeSeq)
#-------------------------------------------
args=options(trailingOnly=TRUE)
data=read.table(args[1])          
#--- raw expression matrix
#
#-------------------------------------------
groups       = rep(c("0hr","4hr","18hr","24hr"),2)
gnames       = rownames(data)
glen         = read.table(args[2]) #--- gene lengths vector
#-------------------------------------------
fit = timeSeq(data, 
              group.label = groups, 
              gene.names  = gnames,
              gene.level  =   TRUE,
              p.values    =   TRUE,
              reads       =   NULL, 
              exon.length =   NULL, 
              gene.length =   glen,
              exon.level  =  FALSE,
              n_cores     =      6,
              iterations  =     10,
              offset      =   TRUE) 
#--- timeseries analysis -------------------
#-------------------------------------------
str(fit)
#-------------------------------------------
NPDE.list =fit$sorted$NPDE_list
NPDE.genes=as.character(NPDE.list[NPDE.list[,2]>=0.6,1])
#-------------------------------------------
PDE.list  =fit$sorted$PDE_list
PDE.genes =as.character(PDE.list[PDE.list[,2]>=0.6,1])
#-------------------------------------------
write.table(NPDE.list,file="./TIMERES_NPDE.txt")
write.table(PDE.list, file="./TIMERES_PDE.txt")
#-------------------------------------------