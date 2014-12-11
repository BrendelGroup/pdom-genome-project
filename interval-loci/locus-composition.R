#!/usr/bin/env Rscript

# Import data
data <- read.table("hymenoptera-ilocus-stats-2014.tsv", header=TRUE, sep="\t")

# Filter by number of genes
data.sub <- data[data$GeneCount == 1,]

# Filter by N content
data.sub <- data.sub[data.sub$Ncontent < 0.25,]

# Filter by length quantiles
qnt <- quantile(data.sub$Length, probs=c(0.1, 0.9))
data.sub <- data.sub[data.sub$Length > qnt[1] & data.sub$Length < qnt[2],]

# Plot Amel vs Pdom
pdom <- data.sub[substr(data.sub$ID, 1, 4) == "Pdom",]
amel <- data.sub[substr(data.sub$ID, 1, 4) == "Amel",]
pdom.h <- hist(pdom$GCcontent, breaks=25, plot=FALSE)
amel.h <- hist(amel$GCcontent, breaks=25, plot=FALSE)
png("amel-pdom-iloci-gc.png", height=1000, width=1000, res=150)
plot(pdom.h$mids, pdom.h$counts, type="l", xlab="%GC content", ylab="Frequency", xlim=c(0,0.65))
lines(amel.h$mids, amel.h$counts, col="orange")
abline(v=median(pdom$GCcontent), lty=3)
abline(v=median(amel$GCcontent), lty=3, col="orange")
d <- dev.off()
