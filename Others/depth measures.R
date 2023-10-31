###Depth measures on the dataset
setwd("C:/Users/marin/OneDrive/Desktop/Erica progetto nps")

dataset<- read.csv("oasis_longitudinal.csv",header=T,sep=',')
dataset <- dataset[-which(is.na(dataset$MMSE)),]
library(MASS)
library(rgl)
library(DepthProc)
library(hexbin)
library(packagefinder)
library(aplpack)
library(robustbase)
#plot confronto tra age e altri parametri + hexagonal binning plot
plot(dataset[,8],dataset[,11], xlab="Age", ylab="MMSE")
bin=hexbin(dataset[,8],dataset[,11], xbins=10, xlab="Age", ylab="MMSE")
plot(bin, main="Hexagonal Binning") #non molto buono visto che i valori di MMSE sono molto concentrati tra 26 e 30

plot(dataset[,8],dataset[,13], xlab="Age", ylab="eTIV")
bin=hexbin(dataset[,8],dataset[,13], xbins=10, xlab="Age", ylab="eTIV")
plot(bin, main="Hexagonal Binning")

#depth measures
dt<-dataset[,c(8,13)]
#turkey depth
tukey_depth=depth(u=dt,method='Tukey')
depth(u = c(77, 1500), X = dt, method = 'Tukey')
depthMedian(dt,depth_params = list(method='Tukey'))
dt[which.max(tukey_depth),]
#max Turkey depth per elemento 112
depthContour(dt,depth_params = list(method='Tukey'))
depthPersp(dt,depth_params = list(method='Tukey'),plot_method = 'rgl')
#Mahalanobis depth
maha_depth <- depth(dt,method='Mahalanobis')
depthMedian(dt,depth_params = list(method='Mahalanobis'))
depthContour(dt,depth_params = list(method='Mahalanobis')) 
depthPersp(dt,depth_params = list(method='Mahalanobis'))
depthPersp(dt,depth_params = list(method='Mahalanobis'),plot_method = 'rgl')

#####per dataset completo
#turkey depth
data<-dataset[,-c(1,2,3,4,6,7)]
data <- data[-which(is.na(data$SES)),]
tukey_depth=depth(u=data,method='Tukey')
depth(u = c(400,77,14,2,28,0.5, 1500,0.7,1), X = data, method = 'Tukey') #provare con valori più indicativi/ di un paziente
depthMedian(data,depth_params = list(method='Tukey'))
data[which.max(tukey_depth),]  
#max Turkey depth per elemento 97

#Mahalanobis depth
maha_depth <- depth(data,method='Mahalanobis')
depthMedian(data,depth_params = list(method='Mahalanobis'))

#outlier detection per dataset bivariato (dt, contenente Age ed eTIV)
depthContour( dt, depth_params = list(method = 'Tukey'), points = TRUE, colors = colorRampPalette(c('white', 'navy')), levels = 10, pdmedian = F, graph_params = list(cex=.01, pch=1), pmean = F )
bagplot(dt)
bagplot_no_outliers <- bagplot(dt) 
outlying_obs <- bagplot_no_outliers$pxy.outlier
#non ci sono outliers 


#outlier detection per dataset completo
bagplot_matrix <- aplpack::bagplot.pairs(data)
#dal plot noto 1/3 outlier. Come li tolgo?
#parte dopo non so come scegliere n_good e n_out nel nostro caso
dt_good <- dt[1:n_good,] 
dt_out <- dt[(n_good+1):n,]
ddPlot(x = dt_good,y = dt_out,depth_params = list(method='Tukey'))