library(mixOmics)

X<-read.csv("polar_pos_pqn_imputed_glog.csv")

sample_data<-read.csv("sample_sheet.csv")

rownames(X)<- X$X

X<-X[,-c(1,2)]

X<-t(X)

#one samples is missing hence  remove the sample from metadata before moving further 

f<-rownames(X)

k<-sample_data$SampleID

missing_sample<-k[!(k %in% f)]

Y<-sample_data[!(sample_data$SampleID %in% missing_sample),]

Y<-Y$REF


#PLS-DA

#plot the samples 

MyResult.plsda <- plsda(X,Y) # without number of variable being selected

plotIndiv(MyResult.plsda)

#customized samples plot

'''This plot depicts the way samples are clustered'''

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

'''Now plot their variable and observe how they cluster'''

plotVar(MyResult.plsda)  

plotVar(MyResult.plsda, var.names=FALSE)

#do a background plot for the samples

background <- background.predict(MyResult.plsda, comp.predicted=2,
                                 dist = "max.dist") 

plotIndiv(MyResult.plsda, comp = 1:2, group = Y,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)









#repeat this for the negative metabolomics data as well. 

z<-read.csv("polar_neg_pqn_imputed_glog.csv")

sample_data<-read.csv("sample_sheet.csv")

rownames(z)<- z$X

z<-z[,-c(1,2)]

z<-t(z)

#remove the missing samples 

f<-rownames(z)

k<-sample_data$SampleID

missing_sample<-k[!(k %in% f)]

Y<-sample_data[!(sample_data$SampleID %in% missing_sample),]

write.csv(Y,"samples_negativemetab.csv")

Y<-Y$REF


#PLS-DA

#plot the samples 

MyResult.plsda <- plsda(z,Y) # without number of variable being selected

plotIndiv(MyResult.plsda)

#customized samples plot

'''This plot depicts the way samples are clustered'''

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

'''Now plot their variable and observe how they cluster'''

plotVar(MyResult.plsda)  

plotVar(MyResult.plsda, var.names=FALSE)

#do a background plot for the samples

plotIndiv(MyResult.plsda, comp = 1:2, group = Y,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)








#repeat the process for RNA seq

rna_seq<-read.csv("rna_norm_counts.csv", header = T, row.names = 1)

rna_seq<-t(rna_seq)

sample_data<-read.csv("sample_sheet.csv")

sample_data<-sample_data$REF



#plot the plsda plot 

MyResult.plsda <- plsda(rna_seq,sample_data) # without number of variable being selected

plotIndiv(MyResult.plsda)


#customized samples plot

'''This plot depicts the way samples are clustered'''

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')



#do a background plot for the samples

background <- background.predict(MyResult.plsda, comp.predicted=2,
                                 dist = "max.dist") 


plotIndiv(MyResult.plsda, comp = 1:2, group = sample_data,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)



'''Now plot their variable and observe how they cluster'''

plotVar(MyResult.plsda)  

plotVar(MyResult.plsda, var.names=FALSE)

'''To many variables are present which makes it confusing to visalize hence
subset the important feature and observe how they cluster'''


#top 10 features selected by RFE method are used to visualize how they cluster 

important_variables<- c("Dapma7bEVm016561","Dapma7bEVm006744","Dapma7bEVm003171","Dapma7bEVm009256",
                        "Dapma7bEVm027566","Dapma7bEVm028337","Dapma7bEVm015665","Dapma7bEVm005195",
                        "Dapma7bEVm023311","Dapma7bEVm012802")

variable_names<-rna_seq[,(colnames(rna_seq) %in% important_variables)]

varible_plot <- plsda(variable_names,sample_data) # without number of variable being selected

plotIndiv(varible_plot)


#plot the variable plot 

plotVar(varible_plot)  

plotVar(varible_plot, var.names=T)


'''This plot depicts the way variables are clustered'''

plotIndiv(varible_plot, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')




















