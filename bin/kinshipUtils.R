#----------------------------------------------------------------------------------
# kinshipUtils
# Michael G. Campana, 2018
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------

bootstrap.kinship <- function(gdsobj, ibdmethod = c("MoM", "MLE"), sample.id = NULL, resample = 100, autosome.only = TRUE, remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, ldmethod = c("composite", "r", "dprime", "corr"), slide.max.bp = 5e+05, slide.max.n = NA, ld.threshold = 0.2, num.thread = 1, verbose = TRUE, kinship = FALSE, kinship.constraint = FALSE, allele.freq = NULL, mlemethod = c("EM", "downhill.simplex", "Jacquard"), max.niter = 1000L, reltol = sqrt(.Machine$double.eps), coeff.correct = TRUE, out.num.iter = TRUE) {
	for (i in 1:resample) {
		snps <- unique(sort(sample(read.gdsn(index.gdsn(gdsobj, "snp.id")), max(read.gdsn(index.gdsn(gdsobj, "snp.id"))), replace=T)))
		snpset <- snpgdsLDpruning(gdsobj, sample.id = sample.id, snp.id = snps, autosome.only = autosome.only, remove.monosnp = remove.monosnp, maf = maf, missing.rate = missing.rate, method = ldmethod, slide.max.bp = slide.max.bp, slide.max.n = slide.max.n, ld.threshold = ld.threshold, verbose = verbose)
		if (ibdmethod == "MoM") {
			snpsetIBD <- snpgdsIBDMoM(gdsobj, sample.id = sample.id, snp.id = unlist(snpset), autosome.only = autosome.only, remove.monosnp = remove.monosnp, maf = maf, missing.rate = missing.rate, kinship = kinship, kinship.constraint = kinship.constraint, num.thread = num.thread, verbose = verbose)
			}
		else {
			snpsetIBD <- snpgdsIBDMLE(gdsobj, sample.id = sample.id, snp.id = unlist(snpset), autosome.only = autosome.only, remove.monosnp = remove.monosnp, maf = maf, missing.rate = missing.rate, kinship = kinship, kinship.constraint = kinship.constraint, allele.freq = allele.freq, method = mlemethod, coeff.correct = coeff.correct, out.num.iter = out.num.iter, num.thread = num.thread, verbose = verbose)
		}
		coeff <- snpgdsIBDSelection(snpsetIBD)
		if (i == 1) {
			res <- cbind(coeff[,1:2],coeff$kinship)
			colnames(res)[i+2] <- paste("Kinship",i, sep="")
			}
		else {
			res <- cbind(res, coeff$kinship)
			colnames(res)[i+2] <- paste("Kinship",i, sep="")
		}
	}
	res <- cbind(res, rowMeans(res[,3:resample+2]))
	sterrs <- apply(res[3:(i+2)], 1, sd)
	sterrs <- sterrs/sqrt(resample)
	res <- cbind(res, sterrs)
	res <- cbind(res, res[i+3] - 1.96 * res[i+4], res[i+3] + 1.96 * res[i+4])
	colnames(res)[i+3] <- "MeanKinship"
	colnames(res)[i+4] <- "StandardError"
	colnames(res)[i+5] <- "LowCI"
	colnames(res)[i+6] <- "HighCI"
	return(res)
}
write.kinship.matrix <- function(x, meanfile = "", cifile = "", sep = ",", digits = 3) {
	n <- length(unique(x$ID1))
	kinship <- format(round(x$MeanKinship, digits = digits), nsmall = digits)
	lowCI <- format(round(x$LowCI, digits = digits), nsmall = digits)
	highCI <- format(round(x$HighCI, digits = digits), nsmall = digits)
	entries <- paste(kinship, " (", lowCI, "-", highCI, ")", sep = "")
	triMatrix <- data.frame(stringsAsFactors = FALSE)
	matStart <- 1
	matEnd <- n
	for (i in 1:n) {
		vals <- c()
		for (j in 1:n) {vals <- c(vals, NA)}
		triMatrix <- rbind(triMatrix, vals)
	}
	rownames(triMatrix) <- unique(x$ID2)
	colnames(triMatrix) <- unique(x$ID1)
	kinMatrix <- triMatrix
	for (i in 1:n) {
		triMatrix[i:n,i] <- entries[matStart:matEnd]
		kinMatrix[i:n,i] <- kinship[matStart:matEnd]
		matStart <- matStart + n - i + 1
		matEnd <- matEnd + n - i
	}
	triMatrix[is.na(triMatrix)] <- ""
	kinMatrix[is.na(kinMatrix)] <- ""
	write.table(triMatrix, file = cifile, sep = sep)
	write.table(kinMatrix, file = meanfile, sep = sep)
}

	
