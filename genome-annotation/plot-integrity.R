data <- read.table("pdom-mrnas-r1.2-qualinfo.tsv",header=TRUE, sep="\t")
png("pdom-genemodel-integrity.png", height=1000, width=1000, res=150)
hist(data$Integrity, breaks=25, border="red", col="#ffeeee", main="",
     xlab="Integrity Score", ylab="# Gene Models")
abline(v=0.95, lwd=2, col="blue", lty=2)
abline(v=0.6, lwd=2, col="blue", lty=2)
d <- dev.off()