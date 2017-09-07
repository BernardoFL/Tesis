library(aire.zmvm)
library(dplyr)
library(reshape2)
library(ggthemes)

o3 <- get_station_data(criterion = "MAXIMOS", 
pollutant = "O3",
year = 2009:2017) 

o3_max <- o3 %>%
group_by(date) %>%
summarise(max = ifelse(all(is.na(value)),
NA,
base::max(value, na.rm = TRUE))) %>%
na.omit()

o3_max$anyo <- sapply(o3_max$date, function(x) substring(x,1,4))
nueve <- o3_max[which(o3_max$anyo == '2009'),]
diez <- o3_max[which(o3_max$anyo == '2010'),]
once <- o3_max[which(o3_max$anyo == '2011'),]
doce <- o3_max[which(o3_max$anyo == '2012'),]
trece <- o3_max[which(o3_max$anyo == '2013'),]
catorce <- o3_max[which(o3_max$anyo == '2014'),]
quince <- o3_max[which(o3_max$anyo == '2015'),]
dieciseis <- o3_max[which(o3_max$anyo == '2016'),]

ggplot(data = nueve) + geom_point(aes(x=1:365,y=max), data=nueve,size=0.03) +
  geom_smooth(aes(x=1:365,y=max), data=diez,size=0.03, alpha=0.1) + geom_smooth(aes(x=1:365,y=max, alpha=0.1), data=once,size=0.03) +
  geom_smooth(aes(x=1:366,y=max, alpha=0.1), data=doce,size=0.03) + geom_smooth(aes(x=1:365,y=max, alpha=0.1), data=trece,size=0.03) +
  geom_smooth(aes(x=1:365,y=max, alpha=0.1), data = catorce,size=0.03) + geom_smooth(aes(x=1:365,y=max, alpha=0.1), data = quince,size=0.03) +
  geom_smooth(aes(x=1:366,y=max, alpha=0.1), data = dieciseis,size=0.03) + theme_bw() + guides(alpha=F)
 # geom_smooth(aes(x=1:366,y=max), data = dieciseis) 

##Prophet

library(prophet)
ds <- c(nueve$date, diez$date, once$date, doce$date, trece$date)
y <- c(nueve$max, diez$max, once$max, doce$max,trece$max)
df <- data.frame(ds,y)
m <- prophet(df)
future <- make_future_dataframe(m,100,freq='day')
pred <- predict(m, future)
plot(m,pred) 


ggplot(data = nueve, aes(x=1:365,y=max)) + geom_point(size=0.3) + 
  geom_smooth(method='gam', formula= y~s(x, bs='ps'))
