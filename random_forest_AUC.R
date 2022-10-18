library(caret)
library(ROSE)
library(DMwR)
library(ggplot2)
library(PRROC)
library(pROC)

###AUC curves are only for a binary classification.
ARGs<-read.csv("C:/Users/usuari/Desktop/usb stick/benjamin/ARGs.csv", sep=";")


dataMatrix <- data.frame(ARGs[-12])

##data partition, for the training dataset
index70 <-createDataPartition(dataMatrix$cultive, p=0.7, list=FALSE)
as.factor(dataMatrix$cultive)
trainSet70 <- dataMatrix[index70,]
testSet70 <- dataMatrix[-index70,]



mtry<-c(1:7)
tunegrid <- expand.grid(mtry=mtry,min.node.size = c(1:7), splitrule="gini")
####for loop to creaate models to mitigate class imbalance (up-sampling, down-sampling, smote, rose) 
modellist <- list()
set.seed(123)
for (i in c("up","down","rose","smote")){
  fitControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             repeats = 3,
                             classProbs =TRUE,
                             ## Evaluate performance with ROC metric using following function (twoClassSummary)
                             summaryFunction = twoClassSummary,
                             preProcOptions = c("center","scale"),
                             search = "grid",
                             sampling=i)

  model_up1<- train(x=trainSet70[,-which(names(trainSet70) %in% c("sample.code","sampletype", "cultive"))],
                    y=trainSet70$cultive, 
                    method="ranger", 
                    metric="ROC", 
                    tuneGrid = tunegrid, 
                    trControl=fitControl, 
                    importance= "permutation")
  key <- toString(i)
  modellist[[key]] <- model_up1
}

###graph comparing the several models accuracy (within each gridline) and between models
results <- resamples(modellist)
summary(modellist)
dotplot(results)

####testmodel
###test models: to get just the AUC,

tests <- list()
for (i in 1:length(modellist)){
predict.upv=predict(modellist[[i]],testSet70[,-which(names(testSet70) %in% c("sample.code","sampletype", "cultive","vegetable"))],type="prob")[,2]
roc.upv <- roc(testSet70$cultive, as.numeric(predict.upv))
rounded_scoresv <- round(as.numeric(predict.upv), digits=1)
roc_roundedupv <- roc(testSet70$cultive, rounded_scoresv)
conf<-auc(roc.upv)
kei <- toString(i)
tests[[kei]] <-conf
}
names(tests)<-names(modellist)

####graph for AUC 
color<-rainbow(length(modellist))
roc_roundedupv<-list()
es<-list()

for (i in 1:length(modellist)){
  predict.upv=predict(modellist[[i]],testSet70[,-which(names(testSet70) %in% c("sample.code","sampletype", "cultive","vegetable"))],type="prob")[,2]
  roc.upv <- roc(testSet70$cultive, as.numeric(predict.upv))
  rounded_scoresv <- round(as.numeric(predict.upv), digits=1)
  roc_roundedupv[[i]] <- roc(testSet70$cultive, rounded_scoresv)
  roc_roundedupv[[i]]$`1-specificity`<-1-roc_roundedupv[[i]]$specificities
  
  #for indivudual plots use the linee below, eliminate lines(), the for loop below and legends
  #plot(y=roc_roundedupv[[i]]$sensitivities, x= roc_roundedupv[[i]]$`1-specificity`, main = "AUC",col="red", type='b',xlab="1-Specificity", ylab="Sensitivity")
  par(pty="s")
  lines(y=roc_roundedupv[[i]]$sensitivities, x= roc_roundedupv[[i]]$`1-specificity`,col=color[i],type='b', lw=1.5)
grid()
  ###list for legend names, exctract the AUC values from list roc_roundedupv generated from previous loop

    for (j in 1:length(roc_roundedupv)){
      es[[j]]<-print(round(roc_roundedupv[[j]]$auc,3))
    }
    legend("bottomright",legend=paste("RF",names(modellist) ,", AUC :",es),bg="white",col=color, lw=2, pt.cex=2, cex=0.73)
    lines(x = c(0,100), y = c(0,100))
    box()
}

title(main="Random Forests AUC",xlab="1-Specificity", ylab="Sensitivity") 
