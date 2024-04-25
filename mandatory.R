library(shiny)
library(ggplot2)
library(tidyverse)
library(data.table)
library(shinydashboard)
library(shinyalert)
library(shinybusy)

library(plyr)
library(reshape2)

PATH = "/Users/jroignant/Desktop/Stage_Martin/complete_app/RNAvis/data/data_reduice/FLEPseq_runs"      #path vers les données réduites
listdir = list.dirs(path = PATH, full.names = FALSE ,recursive = FALSE)                               #list des dossiers dans ce path (list des runs)
PATHo = "/Users/jroignant/Desktop/Stage_Martin/complete_app/RNAvis/data/data/FLEPseq_runs"            #path vers vraies données
barcodeCorrespondancePass = c("")
plotNull = NULL
tableNull = NULL
