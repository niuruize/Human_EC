
################################################################################
# Extended data Fig. 3B
################################################################################
library(randomForestExplainer);library(tidymodels);library(baguette);library(discrim)
# load data
load("~/aging/EC/Human_EC.RData")
# Extracting feature information
DefaultAssay(scRNA) <- "SCT"
sce.markers <- read.csv("markers.txt",sep = "")
Idents(scRNA) <- "Organ"
table(scRNA$Organ)
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(sce.markers$gene)) 
t_expr <- t(as.matrix(sce@assays$SCT@scale.data))
dim(t_expr) 
t_expr[1:4,1:4]
# train data and test data
inTrain <- createDataPartition(y= Idents(sce),p=0.25,list=F)
test_expr <- t_expr[inTrain,]
train_expr <- t_expr[-inTrain,]
test_y <- Idents(sce)[inTrain]
train_y <- Idents(sce)[-inTrain]
save(test_y,train_y, test_expr,train_expr,file = 'input.Rdata') 
train_expr[1:4,1:4]
table(train_y)
table(test_y)
# Building data
test_expr = cbind(data.frame(test_expr),data.frame(Organ=test_y))
train_expr = cbind(data.frame(train_expr),data.frame(Organ=train_y))
rec <- recipe(Organ~.,train_expr)
# 
dt_mod <- decision_tree() %>% set_engine("rpart") %>% set_mode("classification")
kknn_mod <- nearest_neighbor() %>% set_engine('kknn') %>% set_mode('classification')
xgb_mod <- boost_tree() %>% set_engine("xgboost") %>% set_mode("classification")   
nnet_mod <- mlp() %>% set_engine('brulee') %>% set_mode('classification')
naivebayes_mod <- naive_Bayes() %>% set_engine('naivebayes')              
rf_mod <- rand_forest() %>% set_engine('ranger') %>% set_mode('classification')
svm_mod <- svm_rbf() %>% set_engine('kernlab') %>% set_mode('classification')

###
wf <- workflow_set(preproc=list(rec),          
                   models=list(xgb=xgb_mod,          
                               dt=dt_mod,          
                               #log= logistic_mod,          
                               nb=naivebayes_mod,          
                               nnet=nnet_mod,          
                               knn=kknn_mod,          
                               rf=rf_mod,          
                               svm=svm_mod))       
wf

# 
folds <- bootstraps(train_expr,10)
ctr <- control_resamples(save_pred = TRUE)
wf_res <- wf %>% 
  workflow_map("fit_resamples",
               resamples=folds,
               control=ctr)
# results
wf_res
save(wf_res,file = 'wf_res.Rdata') 
#
rank_results(wf_res,rank_metric = "roc_auc") %>% 
  filter(.metric=="roc_auc") %>% 
  select(model,mean)
## visualization
#
pdf("roc_auc.pdf", width = 5,height = 2)
autoplot(wf_res) 

#ggplot2
pdf("auc.pdf", width = 10,height = 10)
collect_predictions(wf_res) %>%           
  group_by(model) %>%  
  roc_curve(Age_Group,.pred_VAT:.pred_Uterus) %>%
  autoplot()
dev.off()

################################################################################
# Extended data Fig. 3C
################################################################################
library(Seurat);library(dplyr);library(randomForest);library(caret);library(pROC);library(caret);library(caret)
DefaultAssay(scRNA) <- "SCT"
sce.markers <- read.csv("markers.txt",sep = "")
Idents(scRNA) <- "Organ"
scRNA <- subset(scRNA, downsample = 2000)
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(sce.markers$gene)) 
t_expr <- t(as.matrix(sce@assays$SCT@scale.data))
dim(t_expr) 
t_expr[1:4,1:4]
table(scRNA$Organ)
# train data and test data
inTrain <- createDataPartition(y= Idents(sce),p=0.25,list=F)
test_expr <- t_expr[inTrain,]
train_expr <- t_expr[-inTrain,]
test_y <- Idents(sce)[inTrain]
train_y <- Idents(sce)[-inTrain]
train_expr[1:4,1:4]
table(train_y)
table(test_y)
save(test_y,train_y, test_expr,train_expr,file = 'input.Rdata') 

#
rf_output = randomForest(x=train_expr, y=train_y,importance = TRUE, ntree = 10001, proximity=TRUE )
rf_output
save(rf_output,file='rf_output.Rdata')

# 75% train data
load(file='rf_output.Rdata')
load(file = 'input.Rdata')
test_outputs <- predict(rf_output,newdata = train_expr,type="prob")
head(test_outputs)
pred_y = colnames(test_outputs)[apply(test_outputs, 1, which.max)]
pred_y = factor(pred_y,levels = levels(train_y))
pdf('RF-performance_train.pdf',width = 8, height = 5)
gplots::balloonplot(table(pred_y,train_y))
dev.off()

# 25% test data
test_outputs <- predict(rf_output,newdata = test_expr,type="prob",importance=TRUE)
save(test_outputs,file='test_output.Rdata')
test_expr[1:4,1:4]
head(test_outputs)
pred_y = colnames(test_outputs)[apply(test_outputs, 1, which.max)]
pred_y = factor(pred_y,levels = levels(test_y))
pdf('RF-performance_test.pdf',width = 10)
gplots::balloonplot(table(pred_y,test_y))
dev.off()
RF.cf <-caret::confusionMatrix(as.factor(pred_y),as.factor(test_y))
RF.cf
save(RF.cf,file='RF_cf.Rdata')

# all data
new_expr <- t(as.matrix(scRNA@assays$RNA@scale.data))
dim(new_expr) 
new_expr[1:4,1:4] 
new_outputs <- predict(rf_output,newdata = new_expr,type="prob")
save(new_outputs,file='rf_output.Rdata')
head(new_outputs)
pred_y = colnames(new_outputs)[apply(new_outputs, 1, which.max)]
new_y <- scRNA$celltype
pred_y = factor(pred_y,levels = levels(new_y))
pdf('RF-performance_all.pdf',width = 8, height = 5)
gplots::balloonplot(table(pred_y,new_y))
dev.off()


################################################################################
# Extended data Fig. 3D
################################################################################
# importance
imp = data.frame(importance(RF.best),MDA.p = RF.best$importanceSD[4])
head(imp)
library(dplyr)
imp = arrange(imp,desc(MeanDecreaseGini)) 
head(imp) 
write.csv(imp,"importance.csv",quote = FALSE)
library(ggplot2)
p1= ggplot(imp,aes(x=MeanDecreaseGini,y=reorder(rownames(imp),MeanDecreaseGini)))+
  geom_bar(position = position_dodge(),
           width = 0.5,
           stat = "identity",
           fill="steelblue")+ 
  theme_minimal() +
  xlab("Mean Decrease in Gini Index")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(size = 16,colour = "black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.x.bottom = element_text(size=16,color="black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #panel.grid = element_blank()
  )
p1
ggsave(filename = "Organ_rf.pdf", plot = p1, device = 'pdf', width = 12, height = 14, units = 'cm')

################################################################################
# Extended data Fig. 3E
################################################################################
load("MacFas_EC.RData")
DefaultAssay(scRNA) <- "RNA"
scRNA = scRNA[,scRNA$Organ %in% c("Adipose","Aorta","Bladder","Breast","Colon","Heart","Kidney","Liver",
                                  "Lung","Muscle","Testis","Trachea","Uterus")]
scRNA$Organ <- as.factor(as.character(scRNA$Organ))
table(scRNA$Organ)
current.cluster.ids <- c("Adipose","Aorta","Bladder","Breast","Colon","Heart","Kidney","Liver","Lung","Muscle","Testis","Trachea","Uterus")
new.cluster.ids <- c("SAT","Aorta","Bladder","MG","Intestines","Heart","Kidney","Liver","Lung","Muscle","Testis","Trachea","Uterus")
scRNA$Organ <- plyr::mapvalues(x = as.character(scRNA$Organ), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA$Organ)
sce.markers <- read.csv("markers.txt",sep = "")
Idents(scRNA) <- "Organ"
scRNA <- subset(scRNA, downsample = 1000)
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(names(rf_output[["forest"]][["ncat"]]))) 
t_expr <- t(as.matrix(sce@assays$RNA@scale.data))
new_expr <- t(as.matrix(sce@assays$RNA@scale.data))
new_expr <- data.frame(new_expr)
names_missing <- names(rf_output[["forest"]][["ncat"]])[!names(rf_output[["forest"]][["ncat"]]) %in% names(new_expr)]
new_expr[,names_missing] <- 0L
new_y <- Idents(sce)
dim(new_expr) 
new_expr[1:4,1:4] 
new_outputs <- predict(rf_output,newdata = as.matrix(new_expr),type="prob")
head(new_outputs)
pred_y = colnames(new_outputs)[apply(new_outputs, 1, which.max)]
pred_y = factor(pred_y,levels = levels(train_y))
pdf('RF-performance_new.pdf',width = 10, height = 6)
gplots::balloonplot(table(pred_y,new_y))+theme(axis.text.x=element_text(vjust = 0, hjust = 0.5))
dev.off()
RF.cf <-caret::confusionMatrix(as.factor(pred_y),as.factor(new_y))
RF.cf

################################################################################
# Extended data Fig. 3F
################################################################################
# homologene
load("Mouse_EC.RData")
DefaultAssay(scRNA) <- "RNA"
library(homologene)
taxData
genelist=rownames(scRNA@assays$RNA@counts)
gene=homologene(genelist, inTax = 10090, outTax = 9606)
gene=gene[!duplicated(gene$`9606`)&!duplicated(gene$`10090`),]
#gene=unique(gene)
table(duplicated(gene$`9606`))
table(duplicated(gene$`10090`))
scRNA <- subset(scRNA, features = gene$`10090`)
rownames(scRNA@assays$RNA@counts)=gene$`9606`
rownames(scRNA@assays$RNA@data)=gene$`9606`
DefaultAssay(scRNA) <- "RNA"
table(scRNA$Organ)
Idents(scRNA) <- "Organ"
scRNA <- subset(scRNA, downsample = 1000)
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(names(rf_output[["forest"]][["ncat"]]))) 
t_expr <- t(as.matrix(sce@assays$RNA@scale.data))
new_expr <- t(as.matrix(sce@assays$RNA@scale.data))
new_expr <- data.frame(new_expr)
names_missing <- names(rf_output[["forest"]][["ncat"]])[!names(rf_output[["forest"]][["ncat"]]) %in% names(new_expr)]
new_expr[,names_missing] <- 0L
new_y <- Idents(sce)
dim(new_expr) 
new_expr[1:4,1:4] 
new_outputs <- predict(rf_output,newdata = as.matrix(new_expr),type="prob")
head(new_outputs)
pred_y = colnames(new_outputs)[apply(new_outputs, 1, which.max)]
pred_y = factor(pred_y,levels = levels(train_y))
pdf('RF-performance_new.pdf',width = 7, height = 10)
gplots::balloonplot(table(new_y,pred_y))+theme(axis.text.x=element_text(vjust = 0, hjust = 0.5))
dev.off()
RF.cf <-caret::confusionMatrix(as.factor(pred_y),as.factor(new_y))
RF.cf








