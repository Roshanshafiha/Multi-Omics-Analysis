
#get the file and edit the data 

X<-read.table(file = 'water_chemicals.tsv', sep = '\t', header = TRUE,row.names = 1)

X<-X[,2:13]

X<-t(X)

site<-1:12

X<-cbind(X,site)

#ensure results are repeatable

set.seed(7)

#load the library

library(mlbench)
library(caret)

# define the control using a random forest selection function

control <- rfeControl(functions=rfFuncs, method="cv", number=10)

# run the RFE algorithm

results <- rfe(X[,1:38],X[,39], sizes=c(1:38), rfeControl=control)

# summarize the results

print(results)

#list the chosen features

top_features<-predictors(results)

# plot the results

plot(results, type=c("g", "o"))




