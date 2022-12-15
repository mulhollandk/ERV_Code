# - - - Full-length HERV

# - Preprocess Samples

herv <- read.csv("hERVsTPM.csv", header=TRUE, row.names=1)
therv <- t(herv)
pheno <- read.csv("Phenotype_tubules.csv", header=TRUE, row.names=1)
rownames(pheno) <- gsub("T","", rownames(pheno))
merge <- merge(therv[,1:2],pheno,by=0)
rownames(merge) <- merge[,1]
pheno <- merge[,c("Fibrosis","Age","Gender","Batch","Race","DM","HTN","RIN","Mito","Dup","Unique","Unmapped")]
phenotype_data_noNA <- pheno[complete.cases(pheno$Fibrosis), ]
merge <- merge(phenotype_data_noNA,therv,by=0)
merge <- as.data.frame(merge)
row.names(merge) <- merge[,1]
merge <- merge[,-1]

# - PCA
pc <- prcomp(merge[, 13: 2610], scale = FALSE)
p1= autoplot(pc, data = merge, colour = "Batch", size = 5, title= "Batch")
p2=autoplot(pc, data = merge, colour = "Gender", size = 5, title= "Gender")
p3=autoplot(pc, data = merge, colour = "Fibrosis", size = 5, title= "Fibrosis")
p4=autoplot(pc, data = merge, colour = "Age", size = 5, title= "Age")
p5=autoplot(pc, data = merge, colour = "Race", size = 5, title= "Race")
p6=autoplot(pc, data = merge, colour = "RIN", size = 5, title= "RIN")
p7=autoplot(pc, data = merge, colour = "Mito", size = 5, title= "Mito")
p8=autoplot(pc, data = merge, colour = "Unmapped", size = 5, title= "Unmapped")
p9=autoplot(pc, data = merge, colour = "Unique", size = 5, title= "Unique")
p10=autoplot(pc, data = merge, colour = "Dup", size = 5, title= "Dup")

# - Linear Regression

GeneName <- colnames(merge[,13: 2610])
GeneName <- as.data.frame(GeneName)
n <- length(merge)
res <- matrix(,n,4)
for (i in c(13: 2610)){
  fit <- lm(merge[,i]~ log2(merge$Fibrosis + 1) + merge$Gender+ merge$Age+ merge$Race+ merge$DM+ merge$HTN+ merge$Batch + merge$RIN + merge$Dup + merge$Mito + merge$Unmapped + merge$Unique)
  res[i,1] <- summary(fit)$coefficients[2,1]
  res[i,2] <- summary(fit)$coefficients[2,4]
}
res[,3] <- p.adjust(res[,2], method="fdr", n=n)
res[,4] <- p.adjust(res[,2], method="bonferroni", n=n)
colnames(res) <- c("Coeff.","Pval","FDR","Bonferroni")
rownames(res) <- colnames(merge)

# - - - Transposable Elements

# - Preprocess Samples

repeat_mask <- read.csv("RepeatMasker_TE_TPM_297.csv", header=TRUE, row.names=1)
targets <- log2(repeat_mask + 1)
#remove outliers
targets  <- subset(targets, select=-c(HK999,HK1943, HK1779, HK1781,HK1783, HK1951, HK1940, HK2244, HK1778, HK1785, HK1784, HK2249, HK1776, HK1833, HK1948))
pheno <- read.csv("Phenotype_tubules.csv", header=TRUE, row.names=1)
rownames(pheno) <- gsub("T","", rownames(pheno))
merge <- merge(targets[,1:2],pheno,by=0)
rownames(merge) <- merge[,1]
pheno <- merge[,c("Fibrosis","Age","Gender","Batch","Race","DM","HTN","RIN","Mito","Dup","Unique","Unmapped")]
phenotype_data_noNA <- pheno[complete.cases(pheno$Fibrosis), ]
merge <- merge(phenotype_data_noNA,therv,by=0)
merge <- as.data.frame(merge)
row.names(merge) <- merge[,1]
merge <- merge[,-1]

# - PCA
pc <- prcomp(merge[, 13: 105880], scale = FALSE)
p1= autoplot(pc, data = merge, colour = "Batch", size = 5, title= "Batch")
p2=autoplot(pc, data = merge, colour = "Gender", size = 5, title= "Gender")
p3=autoplot(pc, data = merge, colour = "Fibrosis", size = 5, title= "Fibrosis")
p4=autoplot(pc, data = merge, colour = "Age", size = 5, title= "Age")
p5=autoplot(pc, data = merge, colour = "Race", size = 5, title= "Race")
p6=autoplot(pc, data = merge, colour = "RIN", size = 5, title= "RIN")
p7=autoplot(pc, data = merge, colour = "Mito", size = 5, title= "Mito")
p8=autoplot(pc, data = merge, colour = "Unmapped", size = 5, title= "Unmapped")
p9=autoplot(pc, data = merge, colour = "Unique", size = 5, title= "Unique")
p10=autoplot(pc, data = merge, colour = "Dup", size = 5, title= "Dup")


# - Linear Regression

GeneName <- colnames(merge[,13: 105880])
GeneName <- as.data.frame(GeneName)
n <- length(merge)
res <- matrix(,n,4)
for (i in c(13: 105880)){
  fit <- lm(merge[,i]~ log2(merge$Fibrosis + 1) + merge$Gender+ merge$Age+ merge$Race+ merge$DM+ merge$HTN+ merge$Batch + merge$RIN + merge$Dup + merge$Mito + merge$Unmapped + merge$Unique)
  res[i,1] <- summary(fit)$coefficients[2,1]
  res[i,2] <- summary(fit)$coefficients[2,4]
}
res[,3] <- p.adjust(res[,2], method="fdr", n=n)
res[,4] <- p.adjust(res[,2], method="bonferroni", n=n)
colnames(res) <- c("Coeff.","Pval","FDR","Bonferroni")
rownames(res) <- colnames(merge)
