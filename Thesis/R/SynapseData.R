
setwd('G:/Thesis Data/OMIC Data/Synapse Data/')
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
install.packages("synapserutils", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapser)
library(synapserutils)
## 
## TERMS OF USE NOTICE:
##   When using Synapse, remember that the terms and conditions of use require that you:
##   1) Attribute data contributors when discussing these data or results from these data.
##   2) Not discriminate, identify, or recontact individuals or groups represented by the data.
##   3) Use and contribute only data de-identified to HIPAA standards.
##   4) Redistribute data only under these same terms of use.
synLogin(email = 'dilarauzuner@gtu.edu.tr',
         password = 'synapsedilara1234')

synLogin('dilarauzuner@gtu.edu.tr', 'synapsedilara1234') 
files <- synapserutils::syncFromSynapse('syn22025006') 
synLogout()
## Welcome, test!
## NULL