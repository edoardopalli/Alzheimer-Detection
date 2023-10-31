#CREAZIONE DATASET

remove(list=ls())
setwd("C:/Users/E5440/Desktop/esami/Nonparametric statistics/Progetto/Alzheimer")

dataset<- read.csv("oasis_longitudinal.csv")
dataset <- dataset[which(dataset$Visit==1),]

#EDUC: years of education 
#SES: socio economic status
#MMSE: mini mental state examination (0:30 30=stato mentale perfetto. Anche persone con un iniziale deterioramento cognitivo, ma con un'alta scolarizzazione possono ottenere un punteggio pari a 29 e 30)
#CDR: clinical dementia rating (0:3 3=perdita di memoria grave)
#eTIV: extimated total intracranial volume 
#nWBV: Normalize Whole Brain Volume
#ASF: Atlas Scaling Factor (the volume-scaling factor required to match each individual to the atlas target)

data <- dataset[,c(6,8,9,11,13,14)]
i.converted <- which(dataset$Group=='Converted')
data <- data[-i.converted,]
data.aux <- dataset[-i.converted,]
data$label <- rep('Nondem', dim(data)[1])
i.dem <- which(data.aux$Group=='Demented')
data$label[i.dem] <- 'Dem'
data$label <- factor(data$label, levels = c('Nondem', 'Dem'))
data$M <- ifelse(data$M.F=='M',1,0)

dataset2 <- read.table('C:/Users/E5440/Downloads/test set.csv', header = TRUE, sep=';')
dataset2$label <- ifelse(dataset2$CDR > 0, 'Dem', 'Nondem')
dataset2  <- dataset2[,c(2,4,5,7,9,10,12)]
dataset2$M <- ifelse(dataset2$M.F=='M',1,0)
colnames(dataset2) <- colnames(data)
data <- rbind(data, dataset2)
write.csv(data, 'data.csv')
write.csv(data[1:271,], 'train.csv')
write.csv(data[272:371,], 'test.csv')

#dataset per regresseione con una variabile
data <- data[,c(1,4,7,8)]
write.csv(data, 'data_unavar.csv')

train <- data[1:271,]
test <- data[272:dim(data)[1],]
write.csv2(test, 'test_unavar.csv')
write.csv2(train, 'train_unavar.csv')


