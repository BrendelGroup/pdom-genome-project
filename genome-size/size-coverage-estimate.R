#!/usr/bin/env Rscript

# Read data into memory
hist17 <- read.table("pdom-17mers.hist", header=FALSE, col.names=c("Matches", "Frequency"))
hist21 <- read.table("pdom-21mers.hist", header=FALSE, col.names=c("Matches", "Frequency"))
hist25 <- read.table("pdom-25mers.hist", header=FALSE, col.names=c("Matches", "Frequency"))
hist29 <- read.table("pdom-29mers.hist", header=FALSE, col.names=c("Matches", "Frequency"))

# Determine the peak frequency in each distribution
max.17 <- max(hist17$Frequency[hist17$Matches > 100 & hist17$Matches < 400])
max.21 <- max(hist21$Frequency[hist21$Matches > 100 & hist21$Matches < 400])
max.25 <- max(hist25$Frequency[hist25$Matches > 100 & hist25$Matches < 400])
max.29 <- max(hist29$Frequency[hist29$Matches > 100 & hist29$Matches < 400])

# Determine the k-mer match value corresponding to the peak frequency in each distribution
# This corresponds to estimated k-mer coverage
cov.17 <- hist17$Matches[hist17$Frequency == max.17]
cov.21 <- hist21$Matches[hist21$Frequency == max.21]
cov.25 <- hist25$Matches[hist25$Frequency == max.25]
cov.29 <- hist29$Matches[hist29$Frequency == max.29]

readcount100bp <- 670323394
readcount35bp  <- 332299160

# Compute preliminary estimates of genome size
ghat.17 <- (readcount100bp * (100 - 17 + 1) + readcount35bp * (35 - 17 + 1)) / cov.17
ghat.21 <- (readcount100bp * (100 - 21 + 1) + readcount35bp * (35 - 21 + 1)) / cov.21
ghat.25 <- (readcount100bp * (100 - 25 + 1) + readcount35bp * (35 - 25 + 1)) / cov.25
ghat.29 <- (readcount100bp * (100 - 29 + 1) + readcount35bp * (35 - 29 + 1)) / cov.29
cat("Preliminary genome size estimates\n")
cat(sprintf("  k=17: %.0f\n", ghat.17))
cat(sprintf("  k=21: %.0f\n", ghat.21))
cat(sprintf("  k=25: %.0f\n", ghat.25))
cat(sprintf("  k=29: %.0f\n", ghat.29))

# Fit linear model of k-mer coverage as a function of k
data <- data.frame(k=c(17,21,25,29), cov=c(cov.17,cov.21,cov.25,cov.29))
mylm <- lm(cov ~ k, data)
slope <- as.numeric(mylm$coefficients[2])
intercept <- as.numeric(mylm$coefficients[1])

# Set k=1, determine coverage and genome size
cov.1 <- slope + intercept
cat(sprintf("Coverage estimate (k=1): %.1f\n", cov.1))
genome.size <- ((readcount100bp * 100)+(readcount35bp * 35)) / cov.1
cat(sprintf("Final genome size estimate (k=1): %.0f\n", genome.size))
cov.13 <- slope*13 + intercept
cov.9  <- slope*9  + intercept
cov.5  <- slope*5  + intercept

# Plot the k-mer distributions
png("pdom-size-kmers.png", height=1000, width=1000, res=150)
plot(hist17$Matches, hist17$Frequency / 1000000, col="red", type="l", xlim=c(1,600), ylim=c(0, 3.5), main="", xlab="K-mer matches", ylab="Frequency (millions)") 
lines(hist21$Matches, hist21$Frequency / 1000000, col="blue")
lines(hist25$Matches, hist25$Frequency / 1000000, col="#009900")
lines(hist29$Matches, hist29$Frequency / 1000000, col="purple")
text(cov.1, .5, sprintf("Coverage = %.1f", cov.1), pos=4, offset=.2)
segments(cov.1,  -.2, cov.1,  .5, lwd=2, col="black")
segments(cov.5,  -.2, cov.5,  .5, lty=3, col="grey")
segments(cov.9,  -.2, cov.9,  .5, lty=3, col="grey")
segments(cov.13, -.2, cov.13, .5, lty=3, col="grey")
abline(v=cov.17, lty=8, col="red")
abline(v=cov.21, lty=8, col="blue")
abline(v=cov.25, lty=8, col="#009900")
abline(v=cov.29, lty=8, col="purple")
d <- dev.off()

# Plot linear model
png("pdom-ck-model.png", height=1000, width=1000, res=150)
plot(data$k, data$cov, main="", xlab="K", ylab="K-mer coverage", xlim=c(0, 30), ylim=c(180, 350))
abline(intercept, slope, col="red", lty=3)
abline(v=1, col="blue", lty=3)
d <- dev.off()

