#!/usr/bin/env Rscript

#----- Import data
#-----------------

# Bees
am <- read.table("amel.seq.dstbn.csv", sep=",", header=FALSE)
bt <- read.table("bter.seq.dstbn.csv", sep=",", header=FALSE)
mr <- read.table("mrot.seq.dstbn.csv", sep=",", header=FALSE)
# Ants
cf <- read.table("cflo.seq.dstbn.csv", sep=",", header=FALSE)
hs <- read.table("hsal.seq.dstbn.csv", sep=",", header=FALSE)
# Wasps
nv <- read.table("nvit.seq.dstbn.csv", sep=",", header=FALSE)
pd <- read.table("pdom.seq.dstbn.csv", sep=",", header=FALSE)

#----- Calculate GC content
#--------------------------

# Bees
am$GCcontent <- 100 * am$V3 / am$V2
bt$GCcontent <- 100 * bt$V3 / bt$V2
mr$GCcontent <- 100 * mr$V3 / mr$V2
# Ants
cf$GCcontent <- 100 * cf$V3 / cf$V2
hs$GCcontent <- 100 * hs$V3 / hs$V2
# Wasps
nv$GCcontent <- 100 * nv$V3 / nv$V2
pd$GCcontent <- 100 * pd$V3 / pd$V2

#----- Select sequences 1 Mbp in length or greater
#-------------------------------------------------

# Bees
amel <- am[am$V4 >= 1000000,]
bter <- bt[bt$V4 >= 1000000,]
mrot <- mr[mr$V4 >= 1000000,]
# Ants
cflo <- cf[cf$V4 >= 1000000,]
hsal <- hs[hs$V4 >= 1000000,]
# Wasps
nvit <- nv[nv$V4 >= 1000000,]
pdom <- pd[pd$V4 >= 1000000,]

#----- Compute histograms oc GC content
#--------------------------------------

# # Bees
# amel.h <- hist(amel$GCcontent, plot=FALSE)
# bter.h <- hist(bter$GCcontent, plot=FALSE)
# mrot.h <- hist(mrot$GCcontent, plot=FALSE)
# # Ants
# cflo.h <- hist(cflo$GCcontent, plot=FALSE)
# hsal.h <- hist(hsal$GCcontent, plot=FALSE)
# # Wasps
# nvit.h <- hist(nvit$GCcontent, plot=FALSE)
# pdom.h <- hist(pdom$GCcontent, plot=FALSE)

#----- Plot GC content
#---------------------

png("hym-gc-dists.png", height=1000, width=1000, res=150)
plot(0, yaxt="n", ylab="", ylim=c(0,0.8), xlim=c(20, 60), xlab="%GC content",
     bty='n')

rug(amel$GCcontent, col="black",  side=3, ticksize=.075, line=-2)
rug(bter$GCcontent, col="black",  side=3, ticksize=.075, line=-4.5)
rug(mrot$GCcontent, col="black",  side=3, ticksize=.075, line=-7)
rug(cflo$GCcontent, col="red",   side=3, ticksize=.075, line=-10.5)
rug(hsal$GCcontent, col="red",   side=3, ticksize=.075, line=-13)
rug(nvit$GCcontent, col="blue",  side=3, ticksize=.075, line=-16.5)
rug(pdom$GCcontent, col="#009900", side=3, ticksize=.075, line=-20)

text(39, 0.73, "A. mellifera",   col="black",  pos=4)
text(43, 0.64, "B. terrestris",  col="black",  pos=4)
text(46, 0.55, "M. rotundata",   col="black",  pos=4)
text(41, 0.42, "C. floridanus",  col="red",   pos=4)
text(53, 0.33, "H. saltator",    col="red",   pos=4)
text(45, 0.21, "N. vitripennis", col="blue",  pos=4)
text(35, 0.08, "P. dominula",    col="#009900", pos=4)

d <- dev.off()
