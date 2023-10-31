setwd('C:/Users/franc/Desktop/NONPA/PROGETTO/ALZHEIMER_prognonpa/EDOARDO')
torte <-  read.csv2('oasis_longitudinal.csv', sep=',')



pie(table(torte$Group),col=c('red','green','blue'))
library(RColorBrewer)
cols <- brewer.pal(3, "RdPu")
pie(table(dataset$Group),col=cols)
unique(dataset$Group)


# GENERALE
cols<- c('#ef756c','#39ba2d','#6b9eff')
pie(table(dataset$Group),col=cols,main='proportion of Demented, Converted & NonDemented')
table(dataset$Group)

# FEMMINE
cols<- c('#ef756c','#39ba2d','#6b9eff')
pie(table(dataset[dataset$M.F=='F',]$Group),col=cols,main='Female pie plot')
table(dataset$Group)

#MASCHI
cols<- c('#ef756c','#39ba2d','#6b9eff')
pie(table(dataset[dataset$M.F=='M',]$Group),col=cols,main='Male pie plot')
table(dataset$Group)
# GENERALE
cols<- c('#ef756c','#39ba2d','#6b9eff')
pie(table(dataset$Group),col=cols,main='')
table(dataset$Group)

# FEMMINE
cols<- c('#ef756c','#39ba2d','#6b9eff')
pie(table(dataset[dataset$M.F=='F',]$Group),col=cols,main='')
table(dataset$Group)

#MASCHI
cols<- c('#ef756c','#39ba2d','#6b9eff')
pie(table(dataset[dataset$M.F=='M',]$Group),col=cols,main='')
table(dataset$Group)




cols<- c('darkblue','red')
pie(table(dataset$Group,dataset$M.F),col=cols,label=lab)
Fem <- brewer.pal(3, "RdPu")
Male <-brewer.pal(3, "Blues")
cols<- c(Fem[1],Fem[2],Fem[3],Male)

lab<-c('F_converted=10','F_Demented=28','F_Nondemented=50','M_converted=4','M_Demented=36','M_Nondemented=22')



