model_gam=gam(label ~ M + s(nWBV,bs='cr') + s(Age, bs='cr')  + s(eTIV, bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam)
pred <- predict(model_gam, newdata = test, type = 'response')
pred <- as.data.frame(pred)
data.roc <- roc.curve(pred, test)

x11()
ggplot(data.roc, aes(x.roc, y.roc))+
  theme(panel.background = element_rect(fill = 'gray90', colour = 'white'))+
  xlim(c(0,1)) + ylim(c(0,1))+
  geom_line(col='dodgerblue4', lwd=1.9)+
  geom_ribbon(aes(ymin=0, ymax=y.roc), fill = 'seagreen2')+
  geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
  xlab('1 - specificity')+ylab('sensitivity')

roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second')#da aggiungere nel plot







