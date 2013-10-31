source("pheatmap.R")

test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# Draw heatmaps
pheatmap(test)

annotation = data.frame(Var1 = factor(1:10 %% 2 == 0,
  labels = c("Class1", "Class2")), Var2 = 1:10)
annotation$Var1 = factor(annotation$Var1, levels = c("Class1", "Class2", "Class3"))
rownames(annotation) = paste("Test", 1:10, sep = "")

pheatmap(test, annotation = annotation)

row_ann <- data.frame(foo=gl(2,nrow(test)/2),`Bar`=relevel(gl(2,nrow(test)/2),"2"))
rownames(row_ann)<-rownames(test)
pheatmap(test, annotation = annotation, annotation_legend = FALSE, drop_levels = FALSE,row_annotation = row_ann)

#Using cytokine annotations
M<-matrix(rnorm(8*20),ncol=8)
row_annotation<-data.frame(A=gl(4,nrow(M)/4),B=gl(4,nrow(M)/4))
eg<-expand.grid(factor(c(0,1)),factor(c(0,1)),factor(c(0,1)))
colnames(eg)<-c("IFNg","TNFa","IL2")
rownames(eg)<-apply(eg,1,function(x)paste0(x,collapse=""))
rownames(M)<-1:nrow(M)
colnames(M)<-rownames(eg)
cytokine_annotation=eg
pheatmap(M,annotation=annotation,row_annotation=row_annotation,annotation_legend=TRUE,row_annotation_legend=TRUE,cluster_rows=FALSE,cytokine_annotation=cytokine_annotation,cluster_cols=FALSE, show_colnames=FALSE)

## Assertions
# Rownames of cytokine annotations == colnames M
stopifnot( rownames(cytokine_annotation) == colnames(M) )

## ... and this is the only requirement?
