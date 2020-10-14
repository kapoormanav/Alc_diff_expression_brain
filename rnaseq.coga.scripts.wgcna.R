ssh -Y interactive6

module load R

library(edgeR)
library(WGCNA)
library(flashClust)
library(DESeq2)

library(parallel)
cl <- makeCluster(4)
library(doParallel)

> exprs <- as.matrix(read.table("RNA.seq.144.coga.inia.count.updated.header.txt", header=TRUE, sep=" ", row.names=1,as.is=TRUE))
> head (exprs)
                X214 X460 X584 X551 X530 X571 X327 X723 X637 X656 X620 X567
ENSG00000223972    1    1    0    0    1    2    0    3    0    0    2    2
ENSG00000227232  200  131   89   92  181  124   77  110   44  115  281  175
ENSG00000243485    0    0    0    0    0    0    0    0    0    0    0    0
ENSG00000237613    0    0    0    0    0    0    0    0    0    0    0    0
ENSG00000268020    0    0    0    0    0    0    0    0    0    0    0    0
ENSG00000240361    0    0    0    0    0    0    0    0    0    0    0    0

cova = read.csv("AUSTRALIAN.BRAIN.PHENO.txt", header=T)
> head (cova)
   IID RNAsequencedby ALC_CAT CON_CASE Age Gender BrainpH  PMI Frozentissue
1  X59           COGA       D  Control  52 Female    6.59 31.0         Left
2  X94           COGA       D  Control  55   Male    6.50 20.0         Left
3  X97           COGA       D  Control  38   Male    6.26 13.5         Left
4  X98           COGA       D  Control  52 Female    5.75  9.5         Left
5 X105           COGA       D  Control  55   Male    6.90  7.5         Left
6 X146           COGA       D  Control  74   Male    6.25 16.5         Left

#Omitting low count rows
dds <- DESeqDataSetFromMatrix(countData = exprs, cova1, design = ~ 1)
dds <- estimateSizeFactors(dds)
isexpr <- rowSums(fpm(dds)>1) >= 0.5 * ncol(dds)
quantLog <- log2( fpm( dds )[isexpr,] + 1)

ddsMat<- DESeqDataSetFromMatrix(exprs, cova1, design=~RNAsequencedby+RIN+CON_CASE)
cova1$CON_CASE <- relevel(cova1$CON_CASE, "Control")
dds<-DESeq(ddsMat, parallel=TRUE)

#map = match( colnames(dat3),rownames(cova3) )
#cova4 = cova3[map,]

#dds = DESeqDataSetFromMatrix( dat3, cova4, design = ~ batch + condition )
#dds = DESeq(dds)

#Omitting low count rows
dds <- rowSums(fpm(dds)>1) >= 0.7 * ncol(dds)
#dds <- dds[ rowSums(counts(dds)) > 1, ]
datExpr0<- assay(dds)

#Checking for good genes
gsg = goodSamplesGenes(datExpr0, verbose = 5);
gsg$allOK

#Varicance Stabilization  & Batch effect Removal
vsd = getVarianceStabilizedData(dds)
design=model.matrix(~sample.table$Type)
vsd.removed=removeBatchEffect(vsd, batch=sample.table$Batch, design= design )
colnames(vsd.removed) <-  dds$SampleNo
vsd.removed.t <- t(vsd.removed)
datExpr<- vsd.removed.t

#Batch effect Removal >>>>>>
isexpr <- rowSums(cpm(exprs)>1) >= 0.5 * ncol(exprs)
gExpr <- DGEList(counts=exprs[isexpr,])
gExpr <- calcNormFactors(gExpr)
design <- model.matrix( ~ RNAsequencedby + RIN, cova1)
vobjGenes <- voom(gExpr, design )
vsd.removed=removeBatchEffect(vobjGenes, batch=cova1$RNAsequencedby, covariates=cova1$RIN)
vsd.removed.t <- t(vsd.removed)
datExpr<- vsd.removed.t

row.names(datExpr) = colnames(exprs)
row.names(cova1) <- rownames(datExpr)
table(rownames(cova1)==rownames(datExpr))

#Clustering of samples
library("amap")
pdf("Sample.clustering.hcluster.of.Batch.pdf")
plot(hcluster(dist(t(vsd.removed)), method= "spearman", link="average" ), cex=1.5, cex.main=2, main= "Hierarchical Clustering of VSD Normalized Count Data ")
dev.off()


# Remove outlying samples 
remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
datExprOut = datExpr[!remove.samples,]
datTraitsOut = cova1[!remove.samples,]
datTraitsOut = datTraitsOut[c(16,17,18)]
dim(datTraitsOut)
save(datExprOut, datTraitsOut, file="SamplesAndTraits_OutliersRemoved.RData")

#Choose a soft threshold power
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExprOut, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

#Plot scale independence
pdf("power.scale.independence.plots.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.90, col="red")

#Plot Mean Connectivity
pdf("power.plots.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

pdf("power.plots.1.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.90, col="red")
dev.off()

softPower = 15
adjacency = adjacency(datExprOut, power = softPower, type = "signed")
TOM = TOMsimilarity(adjacency, TOMType="signed")
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method="average")
pdf("gene.tree.plots.1.pdf")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
dev.off()

#This sets the minimum number of genes to cluster into a module
minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)

dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExprOut, colors= dynamicColors,softPower = 18)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

#Plot eignegenes plots of modules


pdf("MEtree.plots.1.pdf")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
dev.off()

#Merge realted modules
#We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
merge = mergeCloseModules(datExprOut, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file = "geneDendro-3.with.colors.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

###Correlate traits
#Define number of genes and samples
nGenes = ncol(datExprOut)
nSamples = nrow(datExprOut)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExprOut, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraitsOut, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(", 
                    signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))
#display the corelation values with a heatmap plot
pdf(file = "trait.module.corr.pdf", wi = 6.5, he = 6)
labeledHeatmap(Matrix= moduleTraitCor, 
            xLabels= names(datTraits), 
            yLabels= names(MEs), 
            ySymbols= names(MEs), 
            colorLabels= FALSE, 
            colors= blueWhiteRed(50), 
            textMatrix= textMatrix, 
            setStdMargins= FALSE, 
            cex.text= 0.5, 
            zlim= c(-1,1), 
            main= paste("Module-trait relationships"))
            dev.off()



MEs0 = moduleEigengenes(datExprOut, bwColors)$eigengenes
MEs = orderMEs(MEs0)

geneModuleMembership = as.data.frame(cor(datExprOut, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
Alc = as.data.frame(datTraitsOut$alc)
names(Alc) = "alc"
modNames = substring(names(MEs), 3)
names(geneModuleMembership) = paste("MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExprOut, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")


sizeGrWindow(10,7)
pdf(file = "ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = names(datTraitsOut),
             yLabels = names(MEs),
             ySymbols = names(MEs),
             colorLabels = FALSE,
             colors = greenWhiteRed(50),
             textMatrix = textMatrix,
             setStdMargins = FALSE,
             cex.text = 0.5,
             zlim = c(-1,1),
             main= paste("Module-trait relationships"))
                         dev.off()
                         
                         
module = "black"
column = match(module, modNames)
moduleGenes = bwColors==module



sizeGrWindow(7, 7);
pdf(file = "trait.mod.scatter.pdf", wi = 10, he = 7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("Module Membership in", module, "module"),
                 ylab = "Gene significance for alc",
                 main = paste("Module membership vs. gene significance\n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
                 
dev.off()


module = "lightcyan"
column = match(module, modNames)
moduleGenes = bwColors==module



sizeGrWindow(7, 7);
pdf(file = "trait.mod.scatter.pdf", wi = 10, he = 7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("Module Membership in", module, "module"),
                 ylab = "Gene significance for Alc",
                 main = paste("Module membership vs. gene significance\n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
                 
dev.off()



#########################################
R match the columns by rows

https://stackoverflow.com/questions/11977102/order-data-frame-rows-according-to-a-target-vector-that-specifies-the-desired-or


############################################################################################################################

nSamples = nrow(ME.for.audit)
moduleTraitCor = cor(ME.for.audit, audit.scores, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(", 
                    signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(7.5, 8.5, 3, 3))
#display the corelation values with a heatmap plot
pdf(file = "trait.module.corr.pdf", wi = 4.8, he = 6)
labeledHeatmap(Matrix= moduleTraitCor, 
            xLabels= names(audit.scores), 
            yLabels= names(ME.for.audit), 
            ySymbols= names(ME.for.audit), 
            colorLabels= FALSE, 
            colors= blueWhiteRed(50), 
            textMatrix= textMatrix, 
            setStdMargins= FALSE, 
            cex.text= 0.5, 
            zlim= c(-1,1), 
            main= paste("Module-trait relationships"))
            dev.off()
            
            
> names(ME.for.audit)
 [1] "MEturquoise"      "MEdarkolivegreen" "MEdarkturquoise"  "MEblack"         
 [5] "MEwhite"          "MEbisque4"        "MEdarkred"        "MEplum1"         
 [9] "MEthistle1"       "MEdarkorange2"    "MEskyblue3"       "MEbrown"         
[13] "MEskyblue"        "MEbrown4"         "MEroyalblue"      "MEblue"          
[17] "MEdarkgreen"      "MElightcyan1"     "MEplum2"          "MEyellowgreen"   
[21] "MEred"            "MEdarkorange"     "MEivory"          "MEgrey"          
> genes.lightcyan1 = genes[ mergedColors=="lightcyan1"]
> write.table(genes.lightcyan1, "/sc/orga/scratch/kapoom02/genes.in.module.lightcyan1", row.names=FALSE)

genes.black = genes[ mergedColors=="black"]
> write.table(genes.black, "/sc/orga/scratch/kapoom02/genes.in.module.black")
> write.table(genes.black, "/sc/orga/scratch/kapoom02/genes.in.module.black", row.names=FALSE)
> genes.white = genes[ mergedColors=="white"]
> write.table(genes.black, "/sc/orga/scratch/kapoom02/genes.in.module.black", row.names=FALSE)
> write.table(genes.white, "/sc/orga/scratch/kapoom02/genes.in.module.white", row.names=FALSE)
> names(ME.for.audit)
 [1] "MEturquoise"      "MEdarkolivegreen" "MEdarkturquoise"  "MEblack"         
 [5] "MEwhite"          "MEbisque4"        "MEdarkred"        "MEplum1"         
 [9] "MEthistle1"       "MEdarkorange2"    "MEskyblue3"       "MEbrown"         
[13] "MEskyblue"        "MEbrown4"         "MEroyalblue"      "MEblue"          
[17] "MEdarkgreen"      "MElightcyan1"     "MEplum2"          "MEyellowgreen"   
[21] "MEred"            "MEdarkorange"     "MEivory"          "MEgrey"          
> genes.lightcyan1 = genes[ mergedColors=="lightcyan1"]
> write.table(genes.lightcyan1, "/sc/orga/scratch/kapoom02/genes.in.module.lightcyan1", row.names=FALSE)
> dim(ME.for.audit)
[1] 104  24


Residuals:
      Min        1Q    Median        3Q       Max 
-0.150953 -0.053943  0.003595  0.035414  0.286308 

Coefficients:
                              Estimate Std. Error t value Pr(>|t|)    
(Intercept)                 -0.2114964  0.0454195  -4.657 9.91e-06 ***
audit.scores$Classification  0.0322781  0.0150953   2.138   0.0349 *  
audit.scores$Age             0.0028875  0.0006692   4.315 3.76e-05 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.07528 on 100 degrees of freedom
  (1 observation deleted due to missingness)
Multiple R-squared:  0.1864,	Adjusted R-squared:  0.1701 
F-statistic: 11.46 on 2 and 100 DF,  p-value: 3.315e-05



#Batch effect Removal >>>>>>
isexpr <- rowSums(cpm(exprs.audit.matched)>1) >= 0.5 * ncol(exprs)
gExpr <- DGEList(counts=exprs.audit.matched[isexpr,])
gExpr <- calcNormFactors(gExpr)
design <- model.matrix( ~ RNAsequencedby + RIN + Classification, cova3)

form = ~( Age + Gender + Classification + DSMIV + DSMV + PMI + Brain_pH + COD_categoryBrain.pH + Agonal_phase + Liver_class + AUDIT + alcohol_intake_gmsperday + Total_drinking_yrs + RIN + RNAsequencedby )


fit <- lmFit( vobjGenes, model.matrix(~ RNAsequencedby, cova3))
res <- residuals( fit, vobjGenes)

varPartResid <- fitExtractVarPartModel( res, form, cova3 )

plotStratify( Expression ~ Alc_Classification, GE, main=rownames(res)[i])

head (varPartResid[order(-varPartResid$Classification),])
                         Age      Gender Classification        DSMIV
ENSG00000163644 0.0002233958 0.007413620      0.2823613 5.165194e-06
ENSG00000138741 0.1274582313 0.017249953      0.2723797 2.358742e-03
ENSG00000164946 0.0079945864 0.039723943      0.2183398 3.936871e-03
ENSG00000198408 0.0245937740 0.001778879      0.2146053 1.220797e-02
ENSG00000204950 0.0039416132 0.020072534      0.2096210 3.876335e-03
ENSG00000152818 0.0141068075 0.011098754      0.2047287 9.492018e-04
                        DSMV         PMI     Brain_pH COD_category Agonal_phase
ENSG00000163644 1.897264e-03 0.003060275 3.022674e-03 1.213788e-02 1.760725e-05
ENSG00000138741 1.130779e-03 0.011662379 6.692326e-04 3.194120e-02 4.977852e-03
ENSG00000164946 1.069689e-03 0.007618603 1.962210e-02 2.141846e-02 5.303798e-04
ENSG00000198408 5.371706e-03 0.031125613 3.798759e-06 1.021058e-02 6.628984e-02
ENSG00000204950 9.260494e-05 0.002415068 1.008724e-02 2.031436e-02 8.518375e-03
ENSG00000152818 8.689507e-03 0.003174031 1.642756e-02 4.448069e-05 1.407109e-03
                 Liver_class        AUDIT alcohol_intake_gmsperday
ENSG00000163644 0.0012348423 1.597386e-02              0.019700962
ENSG00000138741 0.0008133834 3.586143e-03              0.013299623
ENSG00000164946 0.0016710659 6.198105e-06              0.018689281
ENSG00000198408 0.0233956966 1.920246e-04              0.010421024
ENSG00000204950 0.0002187848 7.537719e-03              0.057866525
ENSG00000152818 0.0004536345 3.740471e-03              0.004923728
                Total_drinking_yrs         RIN RNAsequencedby Residuals
ENSG00000163644        0.002341127 0.055090234   7.076146e-03 0.5884437
ENSG00000138741        0.010891133 0.007935306   1.635814e-05 0.4936300
ENSG00000164946        0.002045277 0.002544972   3.351935e-04 0.6544535
ENSG00000198408        0.010972918 0.015252264   1.180020e-04 0.5734606
ENSG00000204950        0.001727250 0.008031757   1.557285e-06 0.6456773
ENSG00000152818        0.000698335 0.020148678   3.793129e-03 0.7056159




head (moduleColors)
inModule = (moduleColors=="white")
inModule
modProbes = genes[inModule]
head (modProbes)
modTOM = TOM[inModule, inModule]
head (modTOM)
dimnames(modTOM) = list(modProbes, modProbes)
head (annotation)
vis = exportNetworkToVisANT(modTOM,
  file = paste("VisANTInput-", module, ".txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  probeToGene = data.frame(annotation$Gene.stable.ID, annotation$HGNC.symbol) )
module = "white"
vis = exportNetworkToVisANT(modTOM,
  file = paste("VisANTInput-", module, ".txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  probeToGene = data.frame(annotation$Gene.stable.ID, annotation$HGNC.symbol) )
head (vis)
nTop = 30
IMConn = softConnectivity(datExprOut[, modProbes])
  top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
  file = paste("VisANTInput-", module, "-top30.txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  probeToGene = data.frame(annotation$Gene.stable.ID, annotation$HGNC.symbol) )
history (100)


########Volcano plots ###########################

# Download the data from github (click the "raw" button, save as a text file called "results.txt").
# https://gist.github.com/stephenturner/806e31fce55a8b7175af
res <- read.table("results.txt", header=TRUE)
head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))

########Heatmap ###########################

library(“gplots”)
library(“RColorBrewer”)


vsd = varianceStabilizingTransformation(cdsBlind)
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, “GnBu”))(100)
heatmap.2(exprs(vsd)[select,], col = hmcol, trace=”none”, margin=c(10, 6))

########Heatmap concise###########################
dists = dist( t( exprs(vsd) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cdsBlind), paste(condition))
heatmap.2(mat, trace=”none”, col = rev(hmcol), margin=c(13, 13))







