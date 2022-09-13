# "Laboratorio di R per la biostatistica" - Gadda Alberto, Delmiglio Claudia


library(lubridate)
library(readr)
library(randomForest)
library(VSURF)
library(caret)
library(MASS)
library(Boruta)

set.seed(1)

#funzione per fare le confusion matrix
draw_confusion_matrix <- function(cm, title) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(title, cex.main=2)
  
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'CNV', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Other', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'CNV', cex=1.2, srt=90)
  text(140, 335, 'Other', cex=1.2, srt=90)
  
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  






###############################################################
###############################################################
###############################################################
##########                            #########################
##########  DATA CLEANING & EDA       #########################
##########                            #########################
###############################################################
###############################################################
###############################################################

MIOPI_LABORATORIO <- read_delim("C:/Users/gadda/Downloads/MIOPI_LABORATORIO.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)

data = MIOPI_LABORATORIO

sum(is.na(data$`VASCULAR NETWORK on OCTA 1`))
93/122

data=data[,-c(1,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80)] #la prima colonna sono solo indici, mentre l'ultima è composta unicamente da NA
data=data[-c(124:nrow(data)),]

j=1
new_names=c()
for(i in colnames(data)){
  if(substr(i, 1, 1)=="."){
    new_names = c(new_names, as.data.frame(data[1,j])[1,1])
  }
  else{
    new_names = c(new_names, i)
  }
  
  j=j+1
} #alcune intestazioni di colonna si trovano nella prima riga del file excel, mentre alcune nella seconda: creiamo un vettore con i nomi di colonna corretti

colnames(data) = new_names #assegnamo i nuovi nomi di colonna
data = data[-1, ] #eliminiamo la prima riga, ormai inutile

#ORA ABBIAMO UN DATASET CON I NOMI DI COLONNA CORRETTI

#Missing Values
for(i in colnames(data)){
  if(sum(is.na(data[i])) > 0){
    print(c(i, sum(is.na(data[i]))) )
  }
  
}

#Analisi sulle colonne

table(data$SEX)
data$SEX = toupper(data$SEX)
table(data$SEX)
data$SEX = as.factor(data$SEX)


data$AGE = as.numeric(data$AGE)
hist(data$AGE)
boxplot(data$AGE, horizontal = TRUE,notch = TRUE)

data$REF = as.numeric(data$REF)
hist(data$REF)
boxplot(data$REF, horizontal = TRUE,notch = TRUE)

table(data$EYE)
data$EYE = as.factor(data$EYE)


table(data$PSE)
data$PSE = as.factor(data$PSE)

table(data$'RECENT METAMORPHOPSIAS')
data$'RECENT METAMORPHOPSIAS' = as.factor(data$'RECENT METAMORPHOPSIAS')

table(data$POSITION)
data$POSITION = as.factor(data$POSITION)

table(data$haemo)
data$haemo = as.factor(data$haemo)

table(data$'MACULAR CHRA ATROPHY')
data$'MACULAR CHRA ATROPHY' = as.factor(data$'MACULAR CHRA ATROPHY')

table(data$"PAPILLARY CHORIORETINAL ATROPHY")
data$"PAPILLARY CHORIORETINAL ATROPHY" = as.factor(data$"PAPILLARY CHORIORETINAL ATROPHY")

table(data$"ADJACENT CHRA ATROPHY")
data$"ADJACENT CHRA ATROPHY" = as.factor(data$"ADJACENT CHRA ATROPHY")

table(data$"INSIDE CHRA ATROPHY")
data$"INSIDE CHRA ATROPHY" = as.factor(data$"INSIDE CHRA ATROPHY")

table(data$"FUZZY 1")
data$"FUZZY 1" = as.factor(data$"FUZZY 1")

table(data$"FUZZY 2")
data$"FUZZY 2" = as.factor(data$"FUZZY 2")

table(data$"RETINAL thickening 1")
data$"RETINAL thickening 1" = as.factor(data$"RETINAL thickening 1")

table(data$"RETINAL thickening 2")
data$"RETINAL thickening 2" = as.factor(data$"RETINAL thickening 2")

table(data$"IRF 1")
data$"IRF 1" = as.factor(data$"IRF 1")

table(data$"IRF 2")
data$"IRF 2" = as.factor(data$"IRF 2")

table(data$"SRF 1")
data$"SRF 1" = as.factor(data$"SRF 1")

table(data$"SRF 2")
data$"SRF 2" = as.factor(data$"SRF 2")

table(data$"ERM 1")
data$"ERM 1" = as.factor(data$"ERM 1")

table(data$"ERM 2")
data$"ERM 2" = as.factor(data$"ERM 2")

table(data$"ELM PRESENT adjacent to the lesion 1")
data$"ELM PRESENT adjacent to the lesion 1" = as.factor(data$"ELM PRESENT adjacent to the lesion 1")

table(data$"ELM PRESENT adjacent to the lesion 2")
data$"ELM PRESENT adjacent to the lesion 2" = as.factor(data$"ELM PRESENT adjacent to the lesion 2")

elm=c()
for(i in 1:nrow(data)) {
  if((data$"ELM PRESENT within the lesion 1"[i]=="Y") & (data$"ELM INTERRUPTION within the lesion 1"[i]=="Y")){
    elm=c(elm, "I")
  }
  else if((data$"ELM PRESENT within the lesion 1"[i]=="Y") & (data$"ELM INTERRUPTION within the lesion 1"[i]=="N")){
    elm=c(elm, "NI")
  }
  else if(data$"ELM PRESENT within the lesion 1"[i]=="N"){
    elm=c(elm, "N")
  }
  else if(data$"ELM PRESENT within the lesion 1"[i]=="D"){
    elm=c(elm, "D")
  }
  
}
data$"ELM 1" = as.factor(elm)

elm=c()
for(i in 1:nrow(data)) {
  if((data$"ELM PRESENT within the lesion 2"[i]=="Y") & (data$"ELM INTERRUPTION within the lesion 2"[i]=="Y")){
    elm=c(elm, "I")
  }
  else if((data$"ELM PRESENT within the lesion 2"[i]=="Y") & (data$"ELM INTERRUPTION within the lesion 2"[i]=="N")){
    elm=c(elm, "NI")
  }
  else if(data$"ELM PRESENT within the lesion 2"[i]=="N"){
    elm=c(elm, "N")
  }
  else if(data$"ELM PRESENT within the lesion 2"[i]=="D"){
    elm=c(elm, "D")
  }
  
}
data$"ELM 2" = as.factor(elm)

table(data$"EZ PRESENT adjacent to the lesion 1")
data$"EZ PRESENT adjacent to the lesion 1" = as.factor(data$"EZ PRESENT adjacent to the lesion 1")

table(data$"EZ PRESENT adjacent to the lesion 2")
data$"EZ PRESENT adjacent to the lesion 2" = as.factor(data$"EZ PRESENT adjacent to the lesion 2")

ez=c()
for(i in 1:nrow(data)) {
  if((data$"EZ PRESENT within the lesion 1"[i]=="Y") & (data$"EZ INTERRUPTION within the lesion 1"[i]=="Y")){
    ez=c(ez, "I")
  }
  else if((data$"EZ PRESENT within the lesion 1"[i]=="Y") & (data$"EZ INTERRUPTION within the lesion 1"[i]=="N")){
    ez=c(ez, "NI")
  }
  else if(data$"EZ PRESENT within the lesion 1"[i]=="N"){
    ez=c(ez, "N")
  }
  else if(data$"EZ PRESENT within the lesion 1"[i]=="D"){
    ez=c(ez, "D")
  }
  
}
data$"EZ 1" = as.factor(ez)

ez=c()
for(i in 1:nrow(data)) {
  if((data$"EZ PRESENT within the lesion 2"[i]=="Y") & (data$"EZ INTERRUPTION within the lesion 2"[i]=="Y")){
    ez=c(ez, "I")
  }
  else if((data$"EZ PRESENT within the lesion 2"[i]=="Y") & (data$"EZ INTERRUPTION within the lesion 2"[i]=="N")){
    ez=c(ez, "NI")
  }
  else if(data$"EZ PRESENT within the lesion 2"[i]=="N"){
    ez=c(ez, "N")
  }
  else if(data$"EZ PRESENT within the lesion 2"[i]=="D"){
    ez=c(ez, "D")
  }
  
}
data$"EZ 2" = as.factor(ez)                                           
   
data$`HD 1`<-as.factor(data$`HD 1`)
table(data$`HD 1`)

data$`HD 2`<-as.factor(data$`HD 2`)
table(data$`HD 2`)

hd<-table(data$`HD 1`,data$`HD 2`)

data$`SHADOW 1`<-as.factor(data$`SHADOW 1`)
table(data$`SHADOW 1`)

data$`SHADOW 2`<-as.factor(data$`SHADOW 2`)
table(data$`SHADOW 2`)

data$`PERSISTENT CHOROID on SD-OCT 1`<-as.factor(data$`PERSISTENT CHOROID on SD-OCT 1`)
table(data$`PERSISTENT CHOROID on SD-OCT 1`)

data$`PERSISTENT CHOROID on SD-OCT 2`<-as.factor(data$`PERSISTENT CHOROID on SD-OCT 2`)
table(data$`PERSISTENT CHOROID on SD-OCT 2`)

data$`PERSISTENCE OF RPE layer within the lesion 1`<-as.factor(data$`PERSISTENCE OF RPE layer within the lesion 1`)
table(data$`PERSISTENCE OF RPE layer within the lesion 1`)

data$`PERSISTENCE OF RPE layer within the lesion 2`<-as.factor(data$`PERSISTENCE OF RPE layer within the lesion 2`)
table(data$`PERSISTENCE OF RPE layer within the lesion 2`)

data$`RPE PERSISTENT and FLAT 1`<-as.factor(data$`RPE PERSISTENT and FLAT 1`)
table(data$`RPE PERSISTENT and FLAT 1`)

data$`PERSISTENCE OF RPE layer within the lesion 2`<-as.factor(data$`RPE PERSISTENT and FLAT 2`)
table(data$`RPE PERSISTENT and FLAT 2`)

data$`RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 1`<-as.factor(data$`RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 1`)
table(data$`RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 1`)

data$`RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 2`<-as.factor(data$`RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 2`)
table(data$`RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 2`)

data$`RPE PERSISTENT and ELEVATED 1`<-as.factor(data$`RPE PERSISTENT and ELEVATED 1`)
table(data$`RPE PERSISTENT and ELEVATED 1`)

data$`RPE PERSISTENT and ELEVATED 2`<-as.factor(data$`RPE PERSISTENT and ELEVATED 2`)
table(data$`RPE PERSISTENT and ELEVATED 2`)

data$`RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 1`<-as.factor(data$`RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 1`)
table(data$`RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 1`)

data$`RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 2`<-as.factor(data$`RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 2`)
table(data$`RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 2`)

data$`SOPRA LC...57`<-as.factor(data$`SOPRA LC...57`)
table(data$`SOPRA LC...57`)

data$`SOPRA LC...58`<-as.factor(data$`SOPRA LC...58`)
table(data$`SOPRA LC...58`)

data$`VASCULAR NETWORK on OCTA 1`<-as.factor(data$`VASCULAR NETWORK on OCTA 1`)
table(data$`VASCULAR NETWORK on OCTA 1`)

data$`VASCULAR NETWORK on OCTA 2`<-as.factor(data$`VASCULAR NETWORK on OCTA 2`)
table(data$`VASCULAR NETWORK on OCTA 2`)

data$`SCHISIS`<-as.factor(data$`SCHISIS`)
table(data$`SCHISIS`)

data$`INCOMPLETE MH`<-as.factor(data$`INCOMPLETE MH`)
table(data$`INCOMPLETE MH`)

data$`diagnosis 1`<-as.factor(data$`diagnosis 1`)
table(data$`diagnosis 1`)

data$`diagnosis 2`<-as.factor(data$`diagnosis 2`)
table(data$`diagnosis 2`)

final = c()
for(i in data$`FINAL DIAGNOSIS`){
  if(i=="CNV"){final = c(final, "CNV")}
  else{final = c(final, "other")}
}
final
data$`FINAL DIAGNOSIS`=as.factor(final)

leak_ste1<- factor(rep(0,122), levels = c("L","S","N","D"))
for (i in 1:nrow(data)) {
  if(data$`LEAKAGE 1`[i]== "Y"){
    leak_ste1[i] = "L"
  } 
  else if(data$`STAINING 1`[i]=="Y"){
    leak_ste1[i] = "S"
  }
  else if(data$`LEAKAGE 1`[i]== "N" & data$`STAINING 1`[i]=="N"){
    leak_ste1[i] = "N"
  }
  else if(data$`LEAKAGE 1`[i]== "D" & data$`STAINING 1`[i]=="D"){
    leak_ste1[i] = "D"
  }
  else leak_ste1[i] = "NON_TROVATO"
  
}
data$"LEAK_STE 1" = leak_ste1

leak_ste2<- factor(rep(0,122), levels = c("L","S","N","D"))
for (i in 1:nrow(data)) {
  if(data$`LEAKAGE 2`[i]== "Y"){
    leak_ste2[i] = "L"
  } 
  else if(data$`STAINING 2`[i]=="Y"){
    leak_ste2[i] = "S"
  }
  else if(data$`LEAKAGE 2`[i]== "N" && data$`STAINING 2`[i]=="N"){
    leak_ste2[i] = "N"
  }
  else if(data$`LEAKAGE 2`[i]== "D" && data$`STAINING 2`[i]=="D"){
    leak_ste2[i] = "D"
  }
  else leak_ste2[i] = "NON_TROVATO"
  
}
data$"LEAK_STE 2" = leak_ste2

rpe=c()
for(i in 1:nrow(data)) {
  if(data$"PERSISTENCE OF RPE layer within the lesion 1"[i]=="N"){
    rpe=c(rpe, "RPE_N")
  }
  else if(data$"PERSISTENCE OF RPE layer within the lesion 1"[i]=="D"){
    rpe=c(rpe, "RPE_D")
  }
  else if(data$"RPE PERSISTENT and FLAT 1"[i]=="Y"){
    if(data$"RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 1"[i]=="Y"){
      rpe=c(rpe, "RPE_Y_F_Y")
    }
    else{
      rpe=c(rpe, "RPE_Y_F_N")
    }
  }
  else if(data$"RPE PERSISTENT and ELEVATED 1"[i]=="Y"){
    if(data$"RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 1"[i]=="Y"){
      rpe=c(rpe, "RPE_Y_E_Y")
    }
    else{
      rpe=c(rpe, "RPE_Y_E_N")
    }
  }
  else{
    rpe=c(rpe, "RPE_D")
  }
  
}
data$"RPE 1" = as.factor(rpe)

rpe=c()
for(i in 1:nrow(data)) {
  if(data$"PERSISTENCE OF RPE layer within the lesion 2"[i]=="N"){
    rpe=c(rpe, "RPE_N")
  }
  else if(data$"PERSISTENCE OF RPE layer within the lesion 2"[i]=="D"){
    rpe=c(rpe, "RPE_D")
  }
  else if(data$"RPE PERSISTENT and FLAT 2"[i]=="Y"){
    if(data$"RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 2"[i]=="Y"){
      rpe=c(rpe, "RPE_Y_F_Y")
    }
    else{
      rpe=c(rpe, "RPE_Y_F_N")
    }
  }
  else if(data$"RPE PERSISTENT and ELEVATED 2"[i]=="Y"){
    if(data$"RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 2"[i]=="Y"){
      rpe=c(rpe, "RPE_Y_E_Y")
    }
    else{
      rpe=c(rpe, "RPE_Y_E_N")
    }
  }
  else{
    rpe=c(rpe, "RPE_D")
  }
  
}
data$"RPE 2" = as.factor(rpe)



#Elimino le colonne che sono state aggregate in nuove colonne


columns2drop=c("ELM PRESENT within the lesion 1",                        
               "ELM PRESENT within the lesion 2",                       
               "ELM INTERRUPTION within the lesion 1",                   
               "ELM INTERRUPTION within the lesion 2",
               "EZ PRESENT within the lesion 1",                        
               "EZ PRESENT within the lesion 2",                        
               "EZ INTERRUPTION within the lesion 1",                   
               "EZ INTERRUPTION within the lesion 2",
               "LEAKAGE 1",
               "LEAKAGE 2",
               "STAINING 1",
               "STAINING 2",
               "PERSISTENCE OF RPE layer within the lesion 1",          
               "PERSISTENCE OF RPE layer within the lesion 2",          
               "RPE PERSISTENT and FLAT 1",                             
               "RPE PERSISTENT and FLAT 2",                             
               "RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 1",    
               "RPE PERSISTENT and FLAT, WITH FOCAL INTERRUPTION 2",    
               "RPE PERSISTENT and ELEVATED 1",                         
               "RPE PERSISTENT and ELEVATED 2",                         
               "RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 1",
               "RPE PERSISTENT and ELEVATED, WITH FOCAL INTERRUPTION 2",
               "diagnosis 1",
               "diagnosis 2",
               "VASCULAR NETWORK on OCTA 1",
               "VASCULAR NETWORK on OCTA 2"
               )               
data = data[ , -which(names(data) %in% columns2drop)]


#Converto la data in formato date
date_std = function(vect){
  for(x in vect){ 
    if(nchar(x)==4){
      y=paste("01/01/",x, sep="")
      return(dmy(y))
    }
    else {
      return(dmy(x))
    }
  }
  return(ret)
}


d=rep(dmy("1,1,2000"),122)
for(i in c(1:122)){
 d[i]=date_std(data$EVAL[i])
}
data$EVAL = d


#Tolgo gli spazi e i trattini dai nomi delle colonne 
new=c()
for(i in colnames(data)){
  new=c(new, gsub(" ", "_", i, fixed = TRUE))
}
colnames(data)=new

new=c()
for(i in colnames(data)){
  new=c(new, gsub("-", "_", i, fixed = TRUE))
}
colnames(data)=new





data=data[,-3] #elimino la colonna delle date poichè non utilizzabile 

str(data) 






###############################################################
###############################################################
###############################################################
##########                                       ##############
##########  FEATURE SELECTION & CLASSIFICATION   ##############
##########                                       ##############
###############################################################
###############################################################
###############################################################

#FEATURE SELECTION

boruta <- Boruta(FINAL_DIAGNOSIS ~ ., data = data, doTrace = 2, maxRuns = 500)
print(boruta)
plot(boruta, las = 2, cex.axis = 0.7, xlab="" )

importance_boruta = colnames(data)[which(boruta$finalDecision=="Confirmed")]
importance_boruta

data_fs = data[ , which(names(data) %in% c(importance_boruta, "FINAL_DIAGNOSIS"))]


#CLASSIFICATION

train = sample(1:nrow(data), nrow(data)*0.7)


  #RANDOM FOREST

rf_model=randomForest(FINAL_DIAGNOSIS~.,data=data,subset=train, mtry=7)
yhat_rf = predict(rf_model ,newdata=data[-train ,])

cm_rf <- caret::confusionMatrix(data=as.factor(yhat_rf), reference = as.factor(data$FINAL_DIAGNOSIS)[-train])
draw_confusion_matrix(cm_rf, "RANDOM FOREST - Tutte le colonne")

#proviamo con le sole colonne importanti
rf_model_fs=randomForest(FINAL_DIAGNOSIS~.,data=data_fs, subset=train, mtry=7)
yhat_rf_fs = predict(rf_model_fs ,newdata=data_fs[-train ,])

cm_rf_fs <- caret::confusionMatrix(data=as.factor(yhat_rf_fs), reference = as.factor(data_fs$FINAL_DIAGNOSIS)[-train])
draw_confusion_matrix(cm_rf_fs, "RANDOM FOREST - Feature Selection")

random_forest_res = c(cm_rf_fs$overall['Accuracy'],cm_rf_fs$byClass['Sensitivity'],cm_rf_fs$byClass['Specificity'], cm_rf_fs$byClass['Precision'], cm_rf_fs$byClass['Recall'], cm_rf_fs$byClass['F1'] )

  #MLP

library(neuralnet)

bnk_matrix <- model.matrix(~., 
                           data=data)

col_list <- paste(c(colnames(bnk_matrix[,-c(1,62)])),collapse="+")
col_list <- paste(c("FINAL_DIAGNOSISother~",col_list),collapse="")
f <- formula(col_list)

nmodel <- neuralnet(f,data=bnk_matrix[train, ],hidden=c(4),
                    threshold = 0.01,
                    rep=5,
                    learningrate=0.01,
                    learningrate.limit = NULL,
                    learningrate.factor = NULL,
                    algorithm = "rprop+")

output <- compute(nmodel, bnk_matrix[-train,-c(1,62)],rep=1)

pred = as.factor(as.numeric(output$net.result>0.5))
true_label = as.factor(bnk_matrix[-train,62])

cm_mlp <- caret::confusionMatrix(data = pred, reference = true_label)
draw_confusion_matrix(cm_mlp, "MLP - Tutte le colonne")


#proviamo con le sole colonne importanti

bnk_matrix_fs <- model.matrix(~., 
                           data=data_fs)

col_list <- paste(c(colnames(bnk_matrix_fs[,-c(1,19)])),collapse="+")
col_list <- paste(c("FINAL_DIAGNOSISother~",col_list),collapse="")
f <- formula(col_list)

nmodel_fs <- neuralnet(f,data=bnk_matrix_fs[train, ],hidden=c(4),
                    threshold = 0.01,
                    rep=5,
                    learningrate=0.01,
                    learningrate.limit = NULL,
                    learningrate.factor = NULL,
                    algorithm = "rprop+")

output <- compute(nmodel_fs, bnk_matrix_fs[-train,-c(1,19)],rep=1)

pred = as.factor(as.numeric(output$net.result>0.5))
true_label = as.factor(bnk_matrix_fs[-train,19])

cm_mlp_fs <- caret::confusionMatrix(data = pred, reference = true_label)
draw_confusion_matrix(cm_mlp_fs, "MLP - Feature Selection")


#PERFORMANCE
RF = c(cm_rf_fs$overall['Accuracy'],cm_rf_fs$byClass['Sensitivity'],cm_rf_fs$byClass['Specificity'], cm_rf_fs$byClass['Precision'], cm_rf_fs$byClass['Recall'], cm_rf_fs$byClass['F1'] )
MLP = c(cm_mlp_fs$overall['Accuracy'],cm_mlp_fs$byClass['Sensitivity'],cm_mlp_fs$byClass['Specificity'], cm_mlp_fs$byClass['Precision'], cm_mlp_fs$byClass['Recall'], cm_mlp_fs$byClass['F1'] )
names=c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall", "F1")

df <- data.frame(names, RF, MLP)

plot(df$RF, type="b",lwd=3, pch=16,  cex=1.3, col="red", xaxt="n", ann=FALSE, ylim=c(0.87,1.05))
lines(df$MLP, type="o", lwd=3, pch=16, cex=1.5, col="blue")

box()
axis(1, at=1:6, lab=df$names)
legend(1,1.05,c("RandomForest", "MLP"), col=c("red","blue"), lty=c(1,1), pch=c(NA, NA), lwd=c(2,2))






###############################################################
###############################################################
###############################################################
##########                          ###########################
##########  DATA VISUALIZATIONS     ###########################
##########                          ###########################
###############################################################
###############################################################
###############################################################

#TAB1
colnames(data)
library(table1)
label(data$SEX) <- "Sex"
label(data$AGE) <- "Age"
units(data$AGE) <- "years"
tab1<- table1(~ SEX + AGE  | data$FINAL_DIAGNOSIS, data=data, topclass="Rtable1-zebra", overall = "Total" )

#TAB2
colnames(data)
library(table1)
label(data$PERSISTENT_CHOROID_on_SD_OCT_1) <- "PERSISTENT CHROID ON SD OCT 1"
label(data$PERSISTENT_CHOROID_on_SD_OCT_2) <- "PERSISTENT CHROID ON SD OCT 2"
label(data$EZ_2)<- "EZ 2"

tab2<- table1(~ PERSISTENT_CHOROID_on_SD_OCT_1 + PERSISTENT_CHOROID_on_SD_OCT_2 + 
                EZ_2 | data$FINAL_DIAGNOSIS, data=data, topclass="Rtable1-zebra", overall = "Total" )

#grafico1
library(ggplot2)
attach(data)
library(dplyr)
library(ggpubr) 

grafico1 <- ggplot(data, aes(x = FINAL_DIAGNOSIS, y = AGE))
grafico1 + geom_boxplot()
grafico1 + geom_boxplot(notch = TRUE, fill = "lightgray")+
  stat_summary(fun.y = median, geom = "point",
               shape = 18, size = 2.5, color = "#FC4E07")

grafico1<- ggplot(data, aes(x = factor(1), y = AGE)) +
  geom_boxplot(width = 0.4, fill = "lightgray") +
  geom_jitter(aes(color = SEX, shape = FINAL_DIAGNOSIS),
              width = 0.1, size = 2) +
  scale_color_manual(values = c("#E7B800", "#00AFBB")) +
  labs(x = NULL)+
  labs(colour="Sex") +
  labs(shape="Final Diagnosis")+
  scale_y_continuous("Age", limits = c(0,100) )

grafico1

#grafico2_BIS
table(data$FINAL_DIAGNOSIS)
df <- data %>%
  group_by(SEX, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df,4)

grafico2_BIS <- ggplot(df, aes(x = FINAL_DIAGNOSIS, y = counts)) +
  geom_bar(
    aes(color = SEX, fill = SEX),
    stat = "identity", position = position_dodge(0.8),
    width = 0.7,
    
    
  ) +
  scale_color_manual(values = c("#EFC000FF", "#0073C2FF"))+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  scale_x_discrete("Final Diagnosis")+
  scale_y_continuous("Counts", limits = c(0,60))
grafico2_BIS

#grafico2
table(data$FINAL_DIAGNOSIS)
df2 <- data %>%
  group_by(RECENT_METAMORPHOPSIAS, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df2,4)

ap2=data.frame("D","CNV",0)
colnames(ap5)=colnames(df2)
df2 = rbind(df2,ap2 )


grafico2 <- ggplot(df2, aes(x = RECENT_METAMORPHOPSIAS, y = counts)) +
  geom_bar(
    aes(color = FINAL_DIAGNOSIS, fill = FINAL_DIAGNOSIS),
    stat = "identity", position = position_dodge(0.8),
    width = 0.7,
    
    
  ) +
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF"))+
  scale_x_discrete("Recent Metamorphopsias")+
  scale_y_continuous("Counts", limits = c(0,25))
grafico2

#grafico3
df_g3 <- data %>%
  group_by(FUZZY_1, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df_g3,4)

plot(c(1,2,3,4,5))
grafico3<- ggplot(df_g3, aes(x = FUZZY_1 , y = counts)) +
  geom_bar( 
    aes(color = FINAL_DIAGNOSIS, fill = FINAL_DIAGNOSIS),
    stat = "identity", position = position_dodge(0.98),
    width = 0.9,
    
  ) +
  scale_color_manual(values = c("#32cd32", "#ff8c00"))+
  scale_fill_manual(values = c("#32cd32", "#ff8c00"))+
  scale_x_discrete("Fuzzy 1")+
  scale_y_continuous("Counts",limits = c(0,60))

grafico3$labels$fill <- "Final Diagnosis"
grafico3$labels$colour <- "Final Diagnosis"

grafico3

#grafico4
df_g4 <- data %>%
  group_by(FUZZY_2, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df_g4,4)

grafico4<- ggplot(df_g4, aes(x = FUZZY_2 , y = counts)) +
  geom_bar( 
    aes(color = FINAL_DIAGNOSIS, fill = FINAL_DIAGNOSIS),
    stat = "identity", position = position_dodge(0.98),
    width = 0.9,
    
  ) +
  scale_color_manual(values = c("#32cd32", "#ff8c00"))+
  scale_fill_manual(values = c("#32cd32", "#ff8c00"))+
  scale_x_discrete("Fuzzy 2")+
  scale_y_continuous("Counts",limits = c(0,60))

grafico4$labels$fill <- "Final Diagnosis"
grafico4$labels$colour <- "Final Diagnosis"

grafico4

#grafico6
df_g6 <- data %>%
  group_by(SHADOW_1, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df_g6,4)

grafico6<- ggplot(df_g6, aes(x = SHADOW_1 , y = counts)) +
  geom_bar( 
    aes(color = FINAL_DIAGNOSIS, fill = FINAL_DIAGNOSIS),
    stat = "identity", position = position_dodge(0.98),
    width = 0.9,
    
  ) +
  scale_color_manual(values = c("#32cd32", "#ff8c00"))+
  scale_fill_manual(values = c("#32cd32", "#ff8c00"))+
  scale_x_discrete("Shadow 1")+
  scale_y_continuous("Counts",limits = c(0,60))

grafico6$labels$fill <- "Final Diagnosis"
grafico6$labels$colour <- "Final Diagnosis"

grafico6

#grafico7
df_g7 <- data %>%
  group_by(SHADOW_2, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df_g7,4)

grafico7<- ggplot(df_g7, aes(x = SHADOW_2 , y = counts)) +
  geom_bar( 
    aes(color = FINAL_DIAGNOSIS, fill = FINAL_DIAGNOSIS),
    stat = "identity", position = position_dodge(0.98),
    width = 0.9,
    
  ) +
  scale_color_manual(values = c("#32cd32", "#ff8c00"))+
  scale_fill_manual(values = c("#32cd32", "#ff8c00"))+
  scale_x_discrete("Shadow 2")+
  scale_y_continuous("Counts",limits = c(0,60))+
  labs(color = "Final Diagnosis")

grafico7$labels$fill <- "Final Diagnosis"
grafico7$labels$colour <- "Final Diagnosis"

grafico7

#grafico8
df_g8 <- data %>%
  group_by(LEAK_STE_1, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df_g8,4)

ap=data.frame("N","CNV",0)
colnames(ap)=colnames(df_g8)
df_g8 = rbind(df_g8,ap )

#grafico8
grafico8<- ggplot(df_g8, aes(x = LEAK_STE_1 , y =counts)) +
  geom_bar( 
    aes(color = FINAL_DIAGNOSIS, fill = FINAL_DIAGNOSIS),
    stat = "identity", position = position_dodge(0.98),
    width = 0.9
    
  ) +
  scale_color_manual(values = c("#32cd32", "#ff8c00"))+
  scale_fill_manual(values = c("#32cd32", "#ff8c00"))+
  scale_x_discrete("LEAK STE 1")+
  scale_y_continuous("Counts",limits = c(0,65))+
  labs(color = "Final Diagnosis")

grafico8$labels$fill <- "Final Diagnosis"
grafico8$labels$colour <- "Final Diagnosis"

grafico8

#grafico9
df_g9 <- data %>%
  group_by(LEAK_STE_2, FINAL_DIAGNOSIS) %>%
  summarise(counts = n())
head(df_g9,4)

ap2=data.frame("N","CNV",0)
colnames(ap2)=colnames(df_g9)
df_g9 = rbind(df_g9,ap2 )
ap2=data.frame("D","Other",0)
colnames(ap2)=colnames(df_g9)
df_g9 = rbind(df_g9,ap2 )

grafico9<- ggplot(df_g9, aes(x = LEAK_STE_2 , y = counts)) +
  geom_bar( 
    aes(color = FINAL_DIAGNOSIS, fill = FINAL_DIAGNOSIS),
    stat = "identity", position = position_dodge(0.98),
    width = 0.9,
    
  ) +
  scale_color_manual(values = c("#32cd32", "#ff8c00"))+
  scale_fill_manual("Final Diagnosis", values = c("#32cd32", "#ff8c00"))+
  scale_x_discrete("LEAK STE 2")+
  scale_y_continuous("Counts",limits = c(0,65))+
  labs(color = "Final Diagnosis")

grafico9

torta1<-pie(table(data$PERSISTENT_CHOROID_on_SD_OCT_1), col = c("green","blue","yellow"))
torta2<-pie(table(data$PERSISTENT_CHOROID_on_SD_OCT_2), col = c("green","blue","yellow"))
torta3<-pie(table(data$EZ_2), col = c( "green","#00ffff","blue" ,"#ffa500"))











