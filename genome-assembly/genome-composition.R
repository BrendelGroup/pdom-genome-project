#!/usr/bin/env Rscript

# Import data
am <- read.table("amel.seq.dstbn.csv", sep=",", header=FALSE)
bt <- read.table("bter.seq.dstbn.csv", sep=",", header=FALSE)
hs <- read.table("hsal.seq.dstbn.csv", sep=",", header=FALSE)
mr <- read.table("mrot.seq.dstbn.csv", sep=",", header=FALSE)
nv <- read.table("nvit.seq.dstbn.csv", sep=",", header=FALSE)
pd <- read.table("pdom.seq.dstbn.csv", sep=",", header=FALSE)

# Calculate GC content
am$GCcontent <- 100 * am$V3 / am$V2
bt$GCcontent <- 100 * bt$V3 / bt$V2
hs$GCcontent <- 100 * hs$V3 / hs$V2
mr$GCcontent <- 100 * mr$V3 / mr$V2
nv$GCcontent <- 100 * nv$V3 / nv$V2
pd$GCcontent <- 100 * pd$V3 / pd$V2

# Filter out sequences shorter than 1 Mbp
amel <- am[am$V4 >= 1000000,]
bter <- bt[bt$V4 >= 1000000,]
hsal <- hs[hs$V4 >= 1000000,]
mrot <- mr[mr$V4 >= 1000000,]
nvit <- nv[nv$V4 >= 1000000,]
pdom <- pd[pd$V4 >= 1000000,]

# Compute histograms of GC content
amel.h <- hist(amel$GCcontent, plot=FALSE)
bter.h <- hist(bter$GCcontent, plot=FALSE)
hsal.h <- hist(hsal$GCcontent, plot=FALSE)
mrot.h <- hist(mrot$GCcontent, plot=FALSE)
nvit.h <- hist(nvit$GCcontent, plot=FALSE)
pdom.h <- hist(pdom$GCcontent, plot=FALSE)

png("hym-gc-dists.png", height=1000, width=1000, res=150)
plot(0, yaxt="n", ylab="", ylim=c(0,1), xlim=c(20, 60), xlab="%GC content", bty='n')
rug(amel$GCcontent, col="gold",    side=3, ticksize=.1, line=0)
rug(bter$GCcontent, col="#999900", side=3, ticksize=.1, line=-4)
rug(hsal$GCcontent, col="red",     side=3, ticksize=.1, line=-8)
rug(mrot$GCcontent, col="#009900", side=3, ticksize=.1, line=-12)
rug(nvit$GCcontent, col="blue",    side=3, ticksize=.1, line=-16)
rug(pdom$GCcontent, col="black",   side=3, ticksize=.1, line=-20)
text(40, 0.975, "A. mellifera",   col="gold",    pos=4)
text(44, 0.81,  "B. terrestris",  col="#999900", pos=4)
text(37, 0.63,  "H. saltator",    col="red",     pos=2)
text(47, 0.45,  "M. rotundata",   col="#009900", pos=4)
text(45, 0.27,  "N. vitripennis", col="blue",    pos=4)
text(35, 0.095, "P. dominula",    col="black",   pos=4)
d <- dev.off()
