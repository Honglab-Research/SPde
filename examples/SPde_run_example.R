library(forcats)
library(survival)
library(ggplot2)
library(dplyr)
library(survminer)
library(limma)
library(survival)
library(readxl)

#set working directory 

#Load data
aglient_exp1 = read.table('../path/gbm_exp_aglient1.txt')
aglient_exp2= read.table('../path/gbm_exp_aglient2.txt') #you must download from UCSC Xena browser.
gtex40=read.csv('../path/gbm_gtex_40.txt', sep='\t', check.names=F)
rna_exp = read.table('../path/rna_exp.txt', check.names = F)

exp = read.table('../path/agilent_ngs_onlyT_pre_gbm_t_final_exp.txt', check.names = FALSE) #extract only tumor data from output of MergeAgilentNGS.
rtk = read.csv('../path/RTK_list.csv')

annot = readxl::read_excel('../path/Annotation.version2.xlsx') %>% as.data.frame()
pheno = read.table('../path/gbm_pheno.txt', sep='\t', check.names = F)
survival = read.table('../path/TCGA-GBM.survival.tsv.txt', sep='\t', check.names = F)

eff = read.csv('../path/CRISPR_gene_effect.csv')

#mege data
exp <- MergeAgilentNGS(agilent_exp1, agilent_exp2, gtex40, rna_exp)

#Preprocessing
eff <- PrepareCERES(eff)
annot <- PrepareAnnotationFile(annot)

#Run CSGs
csg <- CSGs(eff,annot,exp, pheno, clincal,cancer_type = 'GBM', gef = -1, percent = 1, fc = 1) 

#Run SPde
SP <- SPde(csg$gene, rtk, exp)

###############################################################################################################
# After running this code, you should manually compare the prognosis difference between synergistic pair genes 
# when they act together and when they act individually.
