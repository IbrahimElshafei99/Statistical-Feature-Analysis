#Loading Libraries
library(dplyr)
library(readr) # data input

#Directory
getwd()
setwd("F:/Semester 9/BioStat/Assignment 3/A3")

#Importing data
proteomes <- read.csv("77_cancer_proteomes_CPTAC_itraq.csv", header = TRUE, sep = ";")
clinical <- read.csv("clinical_data_breast_cancer.csv", header = TRUE, sep = ";")
PAM50 <- read.csv("PAM50_proteins.csv", header = TRUE, sep = ";")

#----------------------------------------------------------------------
#Q1
#1
#Remove_NA
proteomes <- proteomes[,-c(2)]
proteomes <- na.omit(proteomes)

#save row_names
row_names <- proteomes$RefSeq_accession_number
row_names

#Transpose necessary col only & drop col(1,2) 
proteomes <- as.data.frame(t(proteomes[,3:85]))
colnames(proteomes) <- row_names


#Append & Rename ID col
proteomes <- cbind(rownames(proteomes), data.frame(proteomes, row.names=NULL))
colnames(proteomes)[1] <- "Complete.TCGA.ID"

#Restructure Function to convert ID format:
convert_ID_format <- function(ID) {                #Input=A2-A0EX.04TCGA
  Part1 = substr(ID, 1, 2)                         #Part1=(A2)
  Part2 = substr(ID, 4, 7)                         #Part2=(A0EX)
  paste("TCGA",Part1,Part2,sep="-")                #Result=TCGA-Part1-Part2="TCGA-A2-A0EX"
}

proteomes$Complete.TCGA.ID <- sapply(proteomes$Complete.TCGA.ID, convert_ID_format)

#Inner join table
New_Data <-  inner_join(clinical, proteomes, by = "Complete.TCGA.ID")

#Remove Duplicates col
which(duplicated(New_Data$Complete.TCGA.ID))
Final_Data <- New_Data [!duplicated(New_Data$Complete.TCGA.ID),]

ncol(clinical)
nrow(Final_Data)


#----------------------------------------------------------------------
#Q2
#1
Final_Data$HER2.Final.Status <- as.integer(as.factor(Final_Data$HER2.Final.Status))
corr <- cor(Final_Data[,31:8024],Final_Data$HER2.Final.Status, method = "pearson")
corr<- as.data.frame(corr)
protein_Names <- rownames(corr)
corr$Protein_ID <- protein_Names

#2
corr <- sapply(corr$V1, abs)
corr <- as.data.frame(corr)
corr$Protein_ID <- protein_Names

order_corr <- corr[order(corr$corr, decreasing = TRUE),]
order_corr <-as.data.frame(order_corr)

#3
summary(order_corr)
High_corr<-order_corr[order_corr$corr>0.4,]
High_corr<-as.data.frame(High_corr)

#--------------------------------------------------------------
#Q3
#1
# Positive_sample<-subset(Final_Data,Final_Data$HER2.Final.Status=="3") 
# Negative_sample<-subset(Final_Data,Final_Data$HER2.Final.Status=="2") 
# Data <- rbind(Positive_sample , Negative_sample)

Data <- subset(Final_Data,Final_Data$HER2.Final.Status!="1") # 1->Equivocal

#2 
t_test_lst = list()

for (i in 31:8024)
{
  t_test <- t.test(Data[,i] ~ HER2.Final.Status , data = Data)
  t_test_lst <- append(t_test_lst, t_test$statistic)
}

t_test_lst = as.data.frame(t(t_test_lst))
t_test_lst$Protein_ID <- protein_Names

#3
t_test_lst <- sapply(t_test_lst$V1, abs)
t_test_lst <- as.data.frame(t_test_lst)
t_test_lst$Protein_ID <- protein_Names

order_t_test_lst <- t_test_lst[order(t_test_lst$t_test_lst, decreasing = TRUE),]
order_t_test_lst <-as.data.frame(order_t_test_lst)

summary(order_t_test_lst)
High_t_test_lst<-order_t_test_lst[order_t_test_lst$t_test_lst>3.5,]
High_t_test_lst<-as.data.frame(High_t_test_lst)

#4
Proteins <-  inner_join(High_corr, High_t_test_lst, by = "Protein_ID")
#High_corr=133 & High_t_test_lst=176 --> common features=78 (almost half)
