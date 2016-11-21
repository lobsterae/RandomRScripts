library(ggplot2)

args=commandArgs(trailingOnly=FALSE)

d   =read.table(args[1])

#---assume input chr - from - to - height

width  = log(1+ifelse(d[,2]>d[,3],d[,2]-d[,3],d[,3]-d[,2]),2)
height = log(1+d[,4],2)

d      = cbind(d,width,height)

pdf("hexprint.pdf")
gplot=ggplot(d, aes(width, height))
gplot + geom_hex()
graphics.off()