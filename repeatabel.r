# --------------------------------------
# Runs a vanilla GWAS in RepeatABEL and plots results.
# --------------------------------------
# Rene Malenfant
# --------------------------------------


# Load libraries
library(RepeatABEL)
library(qqman)
library(qvalue)
#convert.snp.tped(tped = "SLEN/filtered.tped", tfam = "SLEN/filtered.tfam", out = "SLEN/filtered.raw")
#convert.snp.tped(tped = "AXG/filtered.tped", tfam = "AXG/filtered.tfam", out = "AXG/filtered.raw")
#convert.snp.tped(tped = "HdLen/filtered.tped", tfam = "HdLen/filtered.tfam", out = "HdLen/filtered.raw")
#convert.snp.tped(tped = "ZBrd/filtered.tped", tfam = "ZBrd/filtered.tfam", out = "ZBrd/filtered.raw")

# Read & code data
Phen.Data <- read.csv("phenotypes.csv")
Phen.Data$id <- Phen.Data$Bear_ID
Phen.Data$Year <- as.factor(Phen.Data$Year)
Phen.Data$Age2 <- Phen.Data$Age^2
#Phen.Data$Age <- as.factor(Phen.Data$Age)

Phen.Data <- Phen.Data[Phen.Data$Class == "ADULT",]


# SLEN GWAS
# ---------
gen.data <- load.gwaa.data(phe = "dummy_pheno.txt", gen = "SLEN/filtered.raw")
fixed = SLEN ~ Age*Sex
Mod.SLEN <- preFitModel(fixed, random = ~1|id + 1|Year,
                    genabel.data = gen.data, phenotype.data = Phen.Data,
                    corStruc = list(id = list("GRM", "Ind"), Year = list("Ind")))
GWAS.SLEN <- rGLS(fixed, genabel.data = gen.data, phenotype.data = Phen.Data, V = Mod.SLEN$V)

#extract estimated variance components
est.hglm.SLEN <- Mod.SLEN$fitted.hglm
cat("Genotypic, permanent env., and year variance components:\n",
    est.hglm.SLEN$varRanef,", resp.","\n",
    "The residual variance is", est.hglm.SLEN$varFix, ".", "\n")

logVCE <- est.hglm.SLEN$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:\n",
    logVCE[1],"and",logVCE[2],", resp.\n\n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

#print heritability
Vp.SLEN <- est.hglm.SLEN$varFix + est.hglm.SLEN$varRanef[1] + est.hglm.SLEN$varRanef[2] + est.hglm.SLEN$varRanef[3]
cat("h2:", est.hglm.SLEN$varRanef[1]/Vp.SLEN, "\n")
cat("c2:", est.hglm.SLEN$varRanef[2]/Vp.SLEN, "\n")
cat("y2:", est.hglm.SLEN$varRanef[3]/Vp.SLEN, "\n")
cat("e2:", est.hglm.SLEN$varFix/Vp.SLEN, "\n")


# AXG GWAS
# ---------
gen.data <- load.gwaa.data(phe="dummy_pheno.txt", gen="AXG/filtered.raw")
Phen.Data2<-Phen.Data[!is.na(Phen.Data$Num_Cubs),]
fixed = AXG ~ Age*Sex+Age2+D_of_Y+Num_Cubs
Mod.AXG <- preFitModel(fixed, random=~1|id + 1|Year,
                    genabel.data = gen.data, phenotype.data = Phen.Data2,
                    corStruc=list(id=list("GRM","Ind"), Year=list("Ind")))
GWAS.AXG <- rGLS(fixed, genabel.data = gen.data, phenotype.data = Phen.Data2, V = Mod.AXG$V)

#extract estimated variance components
est.hglm.AXG <- Mod.AXG$fitted.hglm
cat("Genotypic, permanent env., and year variance components:","\n",
    est.hglm.AXG$varRanef,", resp.","\n",
    "The residual variance is", est.hglm.AXG$varFix,".","\n")

logVCE <- est.hglm.AXG$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

#print heritability
Vp.AXG <- est.hglm.AXG$varFix + est.hglm.AXG$varRanef[1] + est.hglm.AXG$varRanef[2] + est.hglm.AXG$varRanef[3]
cat("h2:", est.hglm.AXG$varRanef[1]/Vp.AXG, "\n")
cat("c2:", est.hglm.AXG$varRanef[2]/Vp.AXG, "\n")
cat("y2:", est.hglm.AXG$varRanef[3]/Vp.AXG, "\n")
cat("e2:", est.hglm.AXG$varFix/Vp.AXG, "\n")


# HdLen GWAS
# ---------
gen.data <- load.gwaa.data(phe="dummy_pheno.txt", gen="HdLen/filtered.raw")
fixed = HdLen ~ Age*Sex
Mod.HdLen <- preFitModel(fixed, random=~1|id + 1|Year,
                    genabel.data = gen.data, phenotype.data = Phen.Data,
                    corStruc=list(id=list("GRM","Ind"), Year=list("Ind")))
GWAS.HdLen <- rGLS(fixed, genabel.data = gen.data, phenotype.data = Phen.Data, V = Mod.HdLen$V)

#extract estimated variance components
est.hglm.HdLen <- Mod.HdLen$fitted.hglm
cat("Genotypic, permanent env., and year variance components:","\n",
    est.hglm.HdLen$varRanef,", resp.","\n",
    "The residual variance is", est.hglm.HdLen$varFix,".","\n")

logVCE <- est.hglm.HdLen$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

#print heritability
Vp.HdLen <- est.hglm.HdLen$varFix + est.hglm.HdLen$varRanef[1] + est.hglm.HdLen$varRanef[2] + est.hglm.HdLen$varRanef[3]
cat("h2:", est.hglm.HdLen$varRanef[1]/Vp.HdLen, "\n")
cat("c2:", est.hglm.HdLen$varRanef[2]/Vp.HdLen, "\n")
cat("y2:", est.hglm.HdLen$varRanef[3]/Vp.HdLen, "\n")
cat("e2:", est.hglm.HdLen$varFix/Vp.HdLen, "\n")


# ZBrd GWAS
# ---------
gen.data <- load.gwaa.data(phe="dummy_pheno.txt", gen="ZBrd/filtered.raw")
fixed = ZBrd ~ Age*Sex
Mod.ZBrd <- preFitModel(fixed, random=~1|id + 1|Year,
                    genabel.data = gen.data, phenotype.data = Phen.Data,
                    corStruc=list(id=list("GRM","Ind"), Year=list("Ind")))
GWAS.ZBrd <- rGLS(fixed, genabel.data = gen.data, phenotype.data = Phen.Data, V = Mod.ZBrd$V)

#extract estimated variance components
est.hglm.ZBrd <- Mod.ZBrd$fitted.hglm
cat("Genotypic, permanent env., and year variance components:","\n",
    est.hglm.ZBrd$varRanef,", resp.","\n",
    "The residual variance is", est.hglm.ZBrd$varFix,".","\n")

logVCE <- est.hglm.ZBrd$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

#print heritability
Vp.ZBrd <- est.hglm.ZBrd$varFix + est.hglm.ZBrd$varRanef[1] + est.hglm.ZBrd$varRanef[2] + est.hglm.ZBrd$varRanef[3]
cat("h2:", est.hglm.ZBrd$varRanef[1]/Vp.ZBrd, "\n")
cat("c2:", est.hglm.ZBrd$varRanef[2]/Vp.ZBrd, "\n")
cat("y2:", est.hglm.ZBrd$varRanef[3]/Vp.ZBrd, "\n")
cat("e2:", est.hglm.ZBrd$varFix/Vp.ZBrd, "\n")


# Combine & plot

# Choose which SNPs to highlight.  These SNPs were flagged in the paper for being out of HWE, etc.
highlightSNPS <- c("79575_152", "122808_218", "scaffold2_24945886", "scaffold29_6896303", "13214_151", "scaffold3_53165245", "scaffold37_14230470", "scaffold67_3138584", "scaffold93_800986", "112414_185")

GWAS.SLEN.results<-results(GWAS.SLEN)
GWAS.AXG.results<-results(GWAS.AXG)
GWAS.HdLen.results<-results(GWAS.HdLen)
GWAS.ZBrd.results<-results(GWAS.ZBrd)

GWAS.SLEN.results$SNP<-rownames(GWAS.SLEN.results)
GWAS.AXG.results$SNP<-rownames(GWAS.AXG.results)
GWAS.HdLen.results$SNP<-rownames(GWAS.HdLen.results)
GWAS.ZBrd.results$SNP<-rownames(GWAS.ZBrd.results)

GWAS.SLEN.results$Chromosome<-as.integer(levels(GWAS.SLEN.results$Chromosome))[GWAS.SLEN.results$Chromosome]
GWAS.AXG.results$Chromosome<-as.integer(levels(GWAS.AXG.results$Chromosome))[GWAS.AXG.results$Chromosome]
GWAS.HdLen.results$Chromosome<-as.integer(levels(GWAS.HdLen.results$Chromosome))[GWAS.HdLen.results$Chromosome]
GWAS.ZBrd.results$Chromosome<-as.integer(levels(GWAS.ZBrd.results$Chromosome))[GWAS.ZBrd.results$Chromosome]

SetEPS()
postscript("repeatabel.eps")
par(mfrow=c(4,2))
manhattan(GWAS.SLEN.results,chr="Chromosome",bp="Position",p="P1df",highlight=highlightSNPS,genomewideline=-log10(0.05/length(GWAS.SLEN.results$SNP)),suggestiveline=-log10(1/length(GWAS.SLEN.results$SNP)),ylim=c(0,7))
qq(GWAS.SLEN.results$P1df,xlim=c(0,7),ylim=c(0,7))
manhattan(GWAS.AXG.results,chr="Chromosome",bp="Position",p="P1df",highlight=highlightSNPS,genomewideline=-log10(0.05/length(GWAS.AXG.results$SNP)),suggestiveline=-log10(1/length(GWAS.AXG.results$SNP)),ylim=c(0,7))
qq(GWAS.AXG.results$P1df,xlim=c(0,7),ylim=c(0,7))
manhattan(GWAS.HdLen.results,chr="Chromosome",bp="Position",p="P1df",highlight=highlightSNPS,genomewideline=-log10(0.05/length(GWAS.HdLen.results$SNP)),suggestiveline=-log10(1/length(GWAS.HdLen.results$SNP)),ylim=c(0,7))
qq(GWAS.HdLen.results$P1df,xlim=c(0,7),ylim=c(0,7))
manhattan(GWAS.ZBrd.results,chr="Chromosome",bp="Position",p="P1df",highlight=highlightSNPS,genomewideline=-log10(0.05/length(GWAS.ZBrd.results$SNP)),suggestiveline=-log10(1/length(GWAS.ZBrd.results$SNP)),ylim=c(0,7))
qq(GWAS.ZBrd.results$P1df,xlim=c(0,7),ylim=c(0,7))
dev.off()

# Calculate GIFs

estlambda(GWAS.SLEN.results$P1df)$estimate
estlambda(GWAS.AXG.results$P1df)$estimate
estlambda(GWAS.HdLen.results$P1df)$estimate
estlambda(GWAS.ZBrd.results$P1df)$estimate

# Show significant SNPs

GWAS.SLEN@annotation[qvalue(GWAS.SLEN@results$P1df,fdr.level=0.1)$significant,]
GWAS.AXG@annotation[qvalue(GWAS.AXG@results$P1df,fdr.level=0.1)$significant,]
GWAS.HdLen@annotation[qvalue(GWAS.HdLen@results$P1df,fdr.level=0.1)$significant,]
GWAS.ZBrd@annotation[qvalue(GWAS.ZBrd@results$P1df,fdr.level=0.1)$significant,]
