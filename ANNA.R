##--- gene name must be entrez ID 
##-------------------------------
get.ppiNCBI <- function(g.n) 
{
  require(XML)
  ppi <- data.frame()
  for(i in 1:length(g.n)){
    o <- htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/", g.n[i], sep=''))
    # check if interaction table exists
    exist <- length(getNodeSet(o, "//table//th[@id='inter-prod']"))>0
    if(exist){
      p <- getNodeSet(o, "//table")
      ## need to know which table is the good one
      for(j in 1:length(p)){
        int <- readHTMLTable(p[[j]])
        if(colnames(int)[2]=="Interactant"){break}
      }
      ppi <- rbind(ppi, data.frame(egID=g.n[i], intSymbol=int$`Other Gene`))
    }
    # play nice! and avoid being kicked out from NCBI servers
    Sys.sleep(1)
  }
  if(dim(ppi)[1]>0){
    ppi <- unique(ppi)
    print(paste(dim(ppi)[1], "interactions found"))
    return(ppi)
  } else{
    print("No interaction found")
  }
}

library("igraph")
gg <- graph.data.frame(ppi)

plot(gg,
     layout = layout.fruchterman.reingold,
     vertex.label = V(gg)$name,
     vertex.label.color= "black",
     edge.arrow.size=0,
     edge.curved=FALSE
)

tkplot(gg,
## interactive display using tk
       layout = layout.fruchterman.reingold, 
       vertex.label = V(gg)$name,
       vertex.label.color= "black",
       edge.arrow.size=0,
       edge.curved=FALSE
)

#cube roots etc
squidge=function(x,n)
{
  #if n=3 it is a cube root
  d=sign(x)*abs(x)^(1/n)
  d
}
#wss for within group sum squares
wss <- function(d)
{
  sum(scale(d, scale = FALSE)^2)
}

wrap <- function(i, hc, x) 
{
  cl  <- cutree(hc, height=i)
  spl <- split(x, cl)
  wss <- sum(sapply(spl, wss))
  wss
}

# First, import the GTF-file that you have also used as input for htseq-count
library(GenomicFeatures)
# txdb data
txdb = transcriptsBy(Hs.Transcripts,by="gene")
txdb = makeTxDbFromGFF("yourFile.gtf",format="gtf")

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")       
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

transcriptLengths(txdb)

uniformity=function(x)
            {
              d=vector(length=length(x[,1]))
              for(i in 1:length(x[,1]))
                {
                 d[i]=norm(as.matrix(x[i,]),"f")*(sqrt(length(x[i,])-1))/sqrt(length(x[i,])-1) 
                }
              return(d)
             }


sigmoid.fit=nls(~1/(1+((a*(x/0.5)^b)^c)),start=list(a=,b=,c=),x=x,y=y)


#### example simple network plot ####
data

g=graph.edgelist(as.matrix(data[,1:2]),
                 directed=F)

E(g)$weights=data[,3]

plot(g,
     layout=layout.fruchterman.reingold(g,
                                        params=list(weights=E(g)$weights)),
                                        vertex.size=1+2*log(graph.strength(g))*10
     )

# 3d plotting with rgl for data exploration
coords <- layout_with_fr(g, dim=3)
rglplot(g, layout=coords)

#see also daisy for gower dist matrix

hdvm=function(x,y,groups)
{
  x.split=split(x,groups)
  y.split=split(y,groups)
  
  nval=unique(c(x,y))
  sum=0
    
  for(g in 1:length(groups))
  {
    for(i in 1:length(nval))
    {
      nny=length(which(nval[i]==x))
      nnx=length(which(nval[i]==y))
      nx=length(which(nval[i]==x.split[[g]]))
      ny=length(which(nval[i]==y.split[[g]]))
      sum+=abs(nx/nnx-ny/nyy)
    }
  }
  
  return(sum/length(groups))
}

#--- normalized hdvm ----#
nhdvm=function(x,y,groups)
{
  
  d=hdvm(x,y,groups)
  ng=length(groups)
  d=sqrt(ng*pow(d,2))
  return(d)
}

#---fancy curve fitting using the lambda distribution ... a skewed normal

library(GLDEX)
fun.auto.bimodal.ml(data,per.of.mix=0.4,
                    clustering.m=clara,
                    init1.sel="rmfmkl", 
                    init2.sel="rmfmkl",
                    init1=c(0.2,0.5),
                    init2=c(6,0.5),
                    optim.further="Y/N")


GTF2TranscriptDB <- function(gtf.file, out.file = NULL, verbose = TRUE)
{
  require(rtracklayer)
  require(GenomicRanges)
  require(GenomicFeatures)
  
  min.info <- c("gene_id", "transcript_id", "exon_number")
  
  if (verbose) message("Importing ", gtf.file)
  gtf <- import.gff(gtf.file, asRangedData = FALSE)
  
  #parse the per exon attributes
  if (verbose) message("Parsing gene/transcript/exon ids...")
  exon.info <- strsplit(gsub("\"|;", "", values(gtf)$group), split=" ")
  exon.info <- lapply(exon.info, function(x) {
    data <- x[seq(2, length(x), 2)]
    names(data) <- x[seq(1, length(x), 2)]
    data
  })
  
  attribs <- names(exon.info[[1]])
  if (!all(min.info %in% attribs)) stop("Not all required attributes are in this GTF file.")
  exon.info <- do.call(data.frame, c(lapply(attribs, function(x)
    sapply(exon.info, "[", x)),
    stringsAsFactors = FALSE))
  colnames(exon.info) <- attribs
  
  values(gtf) <- exon.info
  
  if (verbose) message("Creating tables...")
  #make transcripts table
  exons.by.tx <- split(gtf, values(gtf)$transcript_id)
  transcripts <- data.frame(
    tx_id = 1:length(exons.by.tx),
    tx_name = names(exons.by.tx),
    tx_chrom = as.character(seqnames(unlist(exons.by.tx))[start(exons.by.tx@partitioning)]),
    tx_strand = as.character(strand(unlist(exons.by.tx))[start(exons.by.tx@partitioning)]),
    tx_start = IRanges::sapply(start(ranges(exons.by.tx)), min),
    tx_end = IRanges::sapply(end(ranges(exons.by.tx)), max),
    stringsAsFactors = FALSE)
  
  #make exons table
  exons.ord <- unlist(exons.by.tx)
  splicings <- data.frame(
    tx_id = rep(1:length(exons.by.tx), elementLengths(exons.by.tx)),
    exon_rank = as.integer(values(exons.ord)$exon_number),
    exon_chrom = as.character(seqnames(exons.ord)),
    exon_strand = as.character(strand(exons.ord)),
    exon_start = start(exons.ord),
    exon_end = end(exons.ord),
    stringsAsFactors = FALSE)
  
  #make genes table
  gene.txs <- tapply(values(gtf)$transcript_id, values(gtf)$gene_id, unique)
  genes <- data.frame(
    tx_name = unlist(gene.txs),
    gene_id = rep(names(gene.txs), sapply(gene.txs, length)),
    stringsAsFactors = FALSE)
  
  #create the db
  if (verbose) message("Creating TranscriptDb. . . . ...")
  gtf.db <- makeTranscriptDb(transcripts, splicings, genes)

  if(!is.null(out.file))
  {
    if (verbose) message("Writing database to file.")
    saveFeatures(gtf.db, out.file)
  }
  gtf.db
}



