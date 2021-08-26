##Gene Co-Expression Analysis 
#Input 
FC_Data <- data.frame(read.csv("celltype1.csv"))
FC_Data1 <- data.frame(read.csv("celltype2.csv"))


##Positively Regulated Cell Type  1 
# Less than 1 log2fold change, but greater than 0 log2fold change
drop1 <- which(FCData$log2FoldChange<1.0 & FCData$log2FoldChange>0 & FCData$padj < 0.05)
within2FC_0FC<- FCData[drop1,] 
dim(within2FC)

# Less than 2 log2fold change, but greater than 1 log2fold change
drop2 <- which(FCData$log2FoldChange<2.0 & FCData$log2FoldChange>1.0 & FCData$padj < 0.05)
within2FC_4FC <- FCData[drop2,] 
dim(within2FC_4FC)
# Less than 3 log2fold change, but greater than 2 log2fold change
drop3 <- which(FCData$log2FoldChange<3 & FCData$log2FoldChange>2.0 & FCData$padj < 0.05)
within4FC_8FC <- FCData[drop3,] 
dim(within4FC_8FC)

##Negatively Regulated Cell Type 1
#Greater than -1 log2fold change, but less than 0 log2fold change
dropNegative1 <- which(FCData$log2FoldChange> -1.0 & FCData$log2FoldChange<0 & FCData$padj < 0.05)
withinNegative2FC_0FC<- FCData[dropNegative1,] 
dim(withinNegative2FC_0FC)

#Greater than -2 log2fold change, but less than -1 log2fold change
dropNegative2 <- which(FCData$log2FoldChange> -2 & FCData$log2FoldChange< -1 & FCData$padj < 0.05)
withinNegative2FC_Negative4FC<- FCData[dropNegative2,] 
dim(withinNegative2FC_Negative4FC)

#Greater than -3 log2fold change, but less than -2 log2fold change
dropNegative3 <- which(FCData$log2FoldChange> -3 & FCData$log2FoldChange< -2 & FCData$padj < 0.05)
withinNegative4FC_Negative8FC<- FCData[dropNegative3,] 
dim(withinNegative4FC_Negative8FC)


##Positively Regulated Cell Type  2
# Less than 1 log2fold change, but greater than 0 log2fold change
drop1 <- which(FCData1$log2FoldChange<1.0 & FCData1$log2FoldChange>0 & FCData1$padj < 0.05)
within2FC_0FC1 <- FCData1[drop1,] 
dim(within2FC1)

# Less than 2 log2fold change, but greater than 1 log2fold change
drop2 <- which(FCData1$log2FoldChange<2.0 & FCData1$log2FoldChange>1.0 & FCData1$padj < 0.05)
within2FC_4FC1 <- FCData1[drop2,] 
dim(within2FC_4FC1)

# Less than 3 log2fold change, but greater than 2 log2fold change
drop3 <- which(FCData1$log2FoldChange<3 & FCData1$log2FoldChange>2.0 & FCData1$padj < 0.05)
within4FC_8FC1 <- FCData1[drop3,] 
dim(within4FC_8FC1)

##Negatively Regulated Cell Type 2
#Greater than -1 log2fold change, but less than 0 log2fold change
dropNegative1 <- which(FCData1$log2FoldChange> -1.0 & FCData1$log2FoldChange<0 & FCData1$padj < 0.05)
withinNegative2FC_0FC1 <- FCData1[dropNegative1,] 
dim(withinNegative2FC_0FC1)

#Greater than -2 log2fold change, but less than -1 log2fold change
dropNegative2 <- which(FCData1$log2FoldChange> -2 & FCData1$log2FoldChange< -1 & FCData1$padj < 0.05)
withinNegative2FC_Negative4FC1 <- FCData1[dropNegative2,] 
dim(withinNegative2FC_Negative4FC1)

#Greater than -3 log2fold change, but less than -2 log2fold change
dropNegative3 <- which(FCData1$log2FoldChange> -3 & FCData1$log2FoldChange< -2 & FCData1$padj < 0.05)
withinNegative4FC_Negative8FC1<- FCData1[dropNegative3,] 
dim(withinNegative4FC_Negative8FC1)


###Merging Gene Expression Data
##Output 
#Co-expressed in 0 to 1 log2fold range 
Co_exp_0to1FC <- merge(within2FC_0FC, within2FC_0FC1, by="Gene")

#Co-expressed in 1 to 2 log2fold range
Co_exp_1to2FC <- merge(within2FC_4FC, within2FC_4FC1, by="Gene")

#Co-expressed in 2 to 3 log2fold range
Co_exp_2to3FC <- merge(within4FC_8FC, within4FC_8FC1, by="Gene")

#Co-expressed in 0 to -1 log2fold range
Co_exp_0toNeg1FC <- merge(withinNegative2FC_0FC, withinNegative2FC_0FC1, by="Gene")
  
#Co-expressed in -1 to -2 log2fold range
Co_exp_Neg1toNeg2FC <- merge(withinNegative2FC_Negative4FC, withinNegative2FC_Negative4FC1, by="Gene")

#Co-expressed in -2 to -3 log2fold range
Co_exp_Neg2toNeg3FC <- merge(withinNegative4FC_Negative8FC, withinNegative4FC_Negative8FC1, by="Gene")


