# --------------------------------------
# Runs a LOCO-style GWAS in RepeatABEL
# For each chromosome or scaffold, a GRM is generated using all SNPs except those
# on the focal chromosome/scaffold. The GWAS for SNPs on that chromosome/scaffold 
# is then conducted. Finally, the results are combined and plotted.
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

GWAS.SLEN.results <- matrix(nrow = 0, ncol = 10)
chrNames <- levels(gen.data@gtdata@chromosome)

for (chr in chrNames) {
  inclSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome == chr]
  inclSNPs.data <- gen.data[,inclSNPs]
  allButSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome != chr]
  allButSNPs.data <- gen.data[,allButSNPs]
  
  allButSNPs.GRM <- compute.GRM(allButSNPs.data)
  
  fixed = SLEN ~ Age*Sex
  Mod.SLEN <- preFitModel(fixed,
                          random = ~1|id + 1|Year,
                          genabel.data = allButSNPs.data,
                          phenotype.data = Phen.Data,
                          GRM = allButSNPs.GRM,
                          corStruc = list(id   = list("GRM", "Ind"),
                                          Year = list("Ind")))
  GWAS.SLEN <- rGLS(fixed,
                    genabel.data = inclSNPs.data,
                    phenotype.data = Phen.Data,
                    GRM = allButSNPs.GRM,
                    V = Mod.SLEN$V)
  GWAS.SLEN.results <- rbind(GWAS.SLEN.results, results(GWAS.SLEN))
}

# AXG GWAS
# ---------
gen.data <- load.gwaa.data(phe="dummy_pheno.txt", gen="AXG/filtered.raw")
Phen.Data2 <- Phen.Data[!is.na(Phen.Data$Num_Cubs),]
GWAS.AXG.results <- matrix(nrow = 0, ncol = 10)
chrNames <- levels(gen.data@gtdata@chromosome)
for (chr in chrNames) {
  inclSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome == chr]
  inclSNPs.data <- gen.data[,inclSNPs]
  allButSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome != chr]
  allButSNPs.data <- gen.data[,allButSNPs]
  
  allButSNPs.GRM <- compute.GRM(allButSNPs.data)
  
  fixed = AXG ~ Age*Sex+Age2+D_of_Y+Num_Cubs
  Mod.AXG <- preFitModel(fixed,
                         random = ~1|id + 1|Year,
                         genabel.data = allButSNPs.data,
                         GRM = allButSNPs.GRM,
                         phenotype.data = Phen.Data2,
                         corStruc = list(id = list("GRM", "Ind"),
                                       Year = list("Ind")))
  GWAS.AXG <- rGLS(fixed,
                   genabel.data = inclSNPs.data,
                   phenotype.data = Phen.Data2,
                   GRM = allButSNPs.GRM,
                   V = Mod.AXG$V)
  GWAS.AXG.results <- rbind(GWAS.AXG.results, results(GWAS.AXG))
}


# HdLen GWAS
# ----------
gen.data <- load.gwaa.data(phe="dummy_pheno.txt", gen="HdLen/filtered.raw")
GWAS.HdLen.results <- matrix(nrow = 0, ncol = 10)
chrNames <- levels(gen.data@gtdata@chromosome)
for (chr in chrNames) {
  inclSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome == chr]
  inclSNPs.data <- gen.data[,inclSNPs]
  allButSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome != chr]
  allButSNPs.data <- gen.data[,allButSNPs]
  
  allButSNPs.GRM <- compute.GRM(allButSNPs.data)
  
  fixed = HdLen ~ Age*Sex
  Mod.HdLen <- preFitModel(fixed,
                           random=~1|id + 1|Year,
                           genabel.data = allButSNPs.data,
                           phenotype.data = Phen.Data,
                           GRM = allButSNPs.GRM,
                           corStruc = list(id = list("GRM","Ind"),
                                           Year = list("Ind")))
  
  GWAS.HdLen <- rGLS(fixed,
                     genabel.data = inclSNPs.data,
                     phenotype.data = Phen.Data,
                     GRM = allButSNPs.GRM,
                     V = Mod.HdLen$V)
  
  GWAS.HdLen.results <- rbind(GWAS.HdLen.results, results(GWAS.HdLen))
}


# ZBrd GWAS
# ----------
gen.data <- load.gwaa.data(phe="dummy_pheno.txt", gen="ZBrd/filtered.raw")
GWAS.ZBrd.results <- matrix(nrow = 0, ncol = 10)
chrNames <- levels(gen.data@gtdata@chromosome)
for (chr in chrNames) {
  inclSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome == chr]
  inclSNPs.data <- gen.data[,inclSNPs]
  allButSNPs <- gen.data@gtdata@snpnames[gen.data@gtdata@chromosome != chr]
  allButSNPs.data <- gen.data[,allButSNPs]
  
  allButSNPs.GRM <- compute.GRM(allButSNPs.data)
  
  fixed = ZBrd ~ Age*Sex
  Mod.ZBrd <- preFitModel(fixed,
                          random=~1|id + 1|Year,
                          genabel.data = allButSNPs.data,
                          phenotype.data = Phen.Data,
                          GRM = allButSNPs.GRM,
                          corStruc = list(id = list("GRM", "Ind"),
                                          Year = list("Ind")))
  
  GWAS.ZBrd <- rGLS(fixed,
                    genabel.data = inclSNPs.data,
                    phenotype.data = Phen.Data,
                    GRM = allButSNPs.GRM,
                    V = Mod.ZBrd$V)
  
  GWAS.ZBrd.results <- rbind(GWAS.ZBrd.results, results(GWAS.ZBrd))
}


# Combine & plot

# Choose which SNPs to highlight.  These SNPs were flagged in the paper for being out of HWE, etc.
highlightSNPS <- c("79575_152", "122808_218", "scaffold2_24945886", "scaffold29_6896303", "13214_151", "scaffold3_53165245", "scaffold37_14230470", "scaffold67_3138584", "scaffold93_800986", "112414_185")

GWAS.SLEN.results$SNP <- rownames(GWAS.SLEN.results)
GWAS.AXG.results$SNP <- rownames(GWAS.AXG.results)
GWAS.HdLen.results$SNP <- rownames(GWAS.HdLen.results)
GWAS.ZBrd.results$SNP <- rownames(GWAS.ZBrd.results)

tmp <- gtdata(gen.data)
GWAS.SLEN.results$Chromosome <- as.integer(as.character(tmp@chromosome[GWAS.SLEN.results$SNP]))
GWAS.AXG.results$Chromosome <- as.integer(as.character(tmp@chromosome[GWAS.AXG.results$SNP]))
GWAS.HdLen.results$Chromosome <- as.integer(as.character(tmp@chromosome[GWAS.HdLen.results$SNP]))
GWAS.ZBrd.results$Chromosome <- as.integer(as.character(tmp@chromosome[GWAS.ZBrd.results$SNP]))

setEPS()
postscript("repeatabel_loco_manhattan.eps")
par(mfrow = c(4, 1))
manhattan(GWAS.SLEN.results, chr="Chromosome", bp="Position", p="P1df",
          highlight = highlightSNPS,
          genomewideline = -log10(0.05/length(GWAS.SLEN.results$SNP)),
          suggestiveline = -log10(1/length(GWAS.SLEN.results$SNP)),
          ylim = c(0, 7))
manhattan(GWAS.AXG.results,
          chr = "Chromosome", bp = "Position", p = "P1df",
          highlight = highlightSNPS,
          genomewideline = -log10(0.05/length(GWAS.AXG.results$SNP)),
          suggestiveline = -log10(1/length(GWAS.AXG.results$SNP)),
          ylim = c(0, 7))
manhattan(GWAS.HdLen.results,
          chr = "Chromosome", bp = "Position", p = "P1df",
          highlight = highlightSNPS,
          genomewideline = -log10(0.05/length(GWAS.HdLen.results$SNP)),
          suggestiveline = -log10(1/length(GWAS.HdLen.results$SNP)),
          ylim = c(0, 7))
manhattan(GWAS.ZBrd.results, chr = "Chromosome", bp = "Position", p = "P1df",
          highlight = highlightSNPS,
          genomewideline = -log10(0.05/length(GWAS.ZBrd.results$SNP)),
          suggestiveline = -log10(1/length(GWAS.ZBrd.results$SNP)),
          ylim = c(0, 7))
dev.off()

setEPS()
postscript("repeatabel_loco_qqplot.eps")
par(mfrow = c(2, 2))
qq(GWAS.SLEN.results$P1df, xlim = c(0, 7), ylim = c(0, 7))
qq(GWAS.AXG.results$P1df, xlim = c(0, 7), ylim = c(0, 7))
qq(GWAS.HdLen.results$P1df, xlim = c(0, 7), ylim = c(0, 7))
qq(GWAS.ZBrd.results$P1df, xlim = c(0, 7),ylim = c(0, 7))
dev.off()

# Calculate GIFs

estlambda(GWAS.SLEN.results$P1df)$estimate
estlambda(GWAS.AXG.results$P1df)$estimate
estlambda(GWAS.HdLen.results$P1df)$estimate
estlambda(GWAS.ZBrd.results$P1df)$estimate

# Show significant SNPs

GWAS.SLEN.results[qvalue(GWAS.SLEN.results$P1df,fdr.level=0.1)$significant,]
GWAS.AXG.results[qvalue(GWAS.AXG.results$P1df,fdr.level=0.1)$significant,]
GWAS.HdLen.results[qvalue(GWAS.HdLen.results$P1df,fdr.level=0.1)$significant,]
GWAS.ZBrd.results[qvalue(GWAS.ZBrd.results$P1df,fdr.level=0.1)$significant,]
