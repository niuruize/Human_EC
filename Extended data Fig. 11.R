
################################################################################
# Extended data Fig. 11C
################################################################################
###----------AutoML-H2O------------###
# Aging prediction models of ECs in different tissues were constructed
library(h2o)
h2o.init()
library(Seurat);library(dplyr);library(randomForest);library(caret);library(pROC);library(caret);library(caret)
library(modeldata);library(tidymodels);library(recipes) 
# load data
load("Brain_Endo.RData")
# Extracting feature information
table(scRNA$Group)
scRNA = scRNA[,scRNA$Group %in% c("Adult","Aging")]
DefaultAssay(scRNA) <- "RNA"
sce.markers <- read.csv("markers.txt",sep = "")
Idents(scRNA) <- "Group"
#sce <- subset(scRNA, subset = nFeature_RNA > 500)
#sce <- subset(sce,downsample = 200)
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(sce.markers$gene)) 
t_expr <- t(as.matrix(sce@assays$RNA@scale.data))
dim(t_expr) 
t_expr[1:4,1:4]

# train data and test data
inTrain <- createDataPartition(y= Idents(sce),p=0.25,list=F)
test_expr <- t_expr[inTrain,]
train_expr <- t_expr[-inTrain,]
test_y <- Idents(sce)[inTrain]
train_y <- Idents(sce)[-inTrain]
#save(test_y,train_y, test_expr,train_expr,file = 'input100.Rdata') 
train_expr[1:4,1:4]
table(train_y)
table(test_y)
table(scRNA$Group)

# Building data
test_expr = cbind(data.frame(test_expr),data.frame(Group=test_y))
train_expr = cbind(data.frame(train_expr),data.frame(Group=train_y))
#train_expr = cbind(data.frame(t_expr),data.frame(Group=Idents(sce)))
rec <- recipe(Group~.,train_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o()
rec$Group <- h2o.asfactor(rec$Group)
y <- "Group"      
x <- setdiff(names(rec),y)

# Building model
am <- h2o.automl(x,y, training_frame = rec,max_models = 20, seed=1234)
save(am,file="h2o_automl.RData")
#am <- h2o.automl(y="Group", training_frame = train_expr, max_models = 10)

#
b <- h2o.get_leaderboard(am)
b
best <- h2o.get_best_model(am)
#best <- h2o.get_best_model(am,"gbm")
perf <- h2o.performance(best)
perf
Confusion_Matrix <- as.matrix(perf@metrics$cm$table[1:2,1:2])
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix1.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')

#
test <- recipe(Group~.,test_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o() 
test$Group <- h2o.asfactor(test$Group)
y <- "Group"
x <- setdiff(names(test),y) 
pred <- h2o.predict(best, test)
##
Confusion_Matrix <- as.matrix(table(cbind(as.data.frame(pred$predict), as.data.frame(test$Group))))
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix2.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')


# ROC-AUC
h2o.auc(perf)
plot(perf,type="roc")
h2o.varimp_plot(best)
#
pdf("permutation_importance.pdf", width = 5,height = 5)
h2o.permutation_importance_plot(best,rec,num_of_features=5)
dev.off()



