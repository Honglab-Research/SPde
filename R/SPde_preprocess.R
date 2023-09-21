#' @title PrepareCERES
#' @description
#' Pre-process CERES data
#'
#' @param eff of CRISPR_gene_effect file
#'
#' @return A data.frame of eff
#'
#' @export
PrepareCERES <- function(eff){
  
  #annotate row and col name
  cells.ind <- eff[[1]]
  rownames(eff) <- cells.ind
  
  col.eff <- colnames(eff)[-1]
  col.eff <- sub("\\.\\.([0-9]+)\\.","_\\1",col.eff)
  
  sym.eff <- sub("_.+","",col.eff)
  eids.eff <- sub(".+_","",col.eff
                  
  eff <- eff[,-1]
  colnames(eff) <- sym.eff 
  eff <- t(eff)

  return(eff)
}

#' @title PrepareAnnotationFile
#' @description
#' Prepare annotation file 
#'
#' @param annot of Annotation.version2.xlsx
#'
#' @return A data.frame of annot
#'
#' @export
#
PrepareAnnotationFile <- function(annot){
  annot = annot[-1,]

  gdsc_depid = annot[which(annot$CRISPR_gene_effect ==1 & annot$GDSC2_public_raw_data_24Jul22==1& annot$`논문_TCGA _annot(from sheet2)` %in%  toupper(cancers_for_eg )& annot$CCLE_expression ==1),]$'DepMap ID'

  ccle_depid = annot[which(annot$CRISPR_gene_effect ==1 & annot$`CCLE_Drug data(CCLE_NP24.2009_Drug_data_2015.02.24)`==1& annot$CCLE_expression ==1 & annot$`논문_TCGA _annot(from sheet2)` %in% toupper(cancers_for_eg )),]$'DepMap ID'

  select = (unique(c(ccle_depid, gdsc_depid)))

  annot = annot[which(annot$`DepMap ID` %in% select),]
  colnames(annot)[1] = 'Depmap_ID' 
  colnames(annot)[22] ='Annotation'
  
  return (annot)
}

#' @title Preprocess
#' @description
#' Pre-process expression data
#'
#' @param exp of merged expression file(TCGA-GBM : NGS, Agilent)
#'
#' @return A data.frame of norm_exp
#'
#' @export
#
preprocess <- function(exp){
  gene = rownames(exp)
  exp = apply(exp,2, as.numeric)
  
  rownames(exp) = gene
  exp = exp[which(rownames(exp)!=""),] 

  #Remove genes for which all samples have a value of 0.
  cnt.zero= rowSums(exp == 0)
  zero_exp = exp[which(cnt.zero != ncol(exp)),]
  
  #Normaliztion
  norm_exp = t(scale(t(zero_exp)))
  
  return(norm_exp)
}

#' @title deg
#' @description
#' DEG anlaysis using limma 
#'
#' @param exp
#' @param metadata of TCGA-GBM phenotype
#'
#' @return list of exp_TN(preprocessed expression data), Fgenes(Final significantly differental expressed genes) , sam(sampleID used in anlysis)
#'
#' @export
#
deg <- function(exp,metadata,fc){
  metadata = metadata[which(metadata$sampleID %in% colnames(exp)),c('sampleID','sample_type')]
  
  sam = metadata[which(metadata$sample_type != ''),]$sampleID
  exp_TN = exp[,sam] 
  
  rownames(metadata) = metadata$sampleID
  metadata = (metadata[colnames(exp_TN),]) 
  
  only_N = metadata[which(metadata$sample_type == 'Solid Tissue Normal'),]$sampleID
  group = sapply(colnames(exp_TN), function(x) ifelse(x %in% only_N, 'Normal','Tumor')) %>% as.data.frame()
  rownames(group) = colnames(exp_TN)  
  
  #limma
  design <- model.matrix(~0+group$.)
  colnames(design) <-c('Normal','Tumor') #check always 
  fit <- lmFit(exp_TN, design)
  cont <- makeContrasts(Tumor-Normal, levels=design)
  fit.cont <- contrasts.fit(fit, cont)
  fit.cont <- eBayes(fit.cont)
  res <- topTable(fit.cont, number = Inf)
  dim(res) 

  #Extract significantly differential expressed genes between normal and tumor
  upres = res[which(res$logFC >= fc),]
  dnres = res[which(res$logFC <= -fc),]
  f_res = rbind(upres,dnres)
  
  g= rownames(f_res[which(f_res$P.Value <=0.05),])
  f_res = f_res[g,]
  
  Fgenes = rownames(f_res)
  exp_TN = as.data.frame(exp_TN)
  colnames(exp_TN) = sam
  
  return(list(exp_TN,Fgenes,sam))
}
                 
#' @title MergeAgilentNGS
#' @description
#' Merge agilent, NGS, GTEx data.                
#' Normalize each data, preprocess : scale 
#'
#' @param agilent1; Xena browser provides TCGA-GBM Agilent data separately.
#' @param agilent1; Xena browser provides TCGA-GBM Agilent data separately.
#' @param gtex40; randomly extracted brain normal data from GTEx db.
#' @param rna_exp; TCGA-GBM expression data(NGS).
#'
#' @return A data frame of test_gbm. 
#'
#' @export
#

MergeAgilentNGS <- function(agilent1, agilent2, gtex40, rna_exp){
  #aglinat1
  colnames(agilent_exp1) = agilent_exp1[1,]
  rownames(agilent_exp1) = agilent_exp1[,1]
  agilent_exp1 = agilent_exp1[-1,-1]

  gene = rownames(agilent_exp1)
  agilent_exp1 = apply(agilent_exp1,2,as.numeric)
  rownames(agilent_exp1) = gene

  #agilent2
  colnames(agilent_exp2) = agilent_exp2[1,]
  rownames(agilent_exp2) = agilent_exp2[,1]
  agilent_exp2 = agilent_exp2[-1,-1]
  gene = rownames(agilent_exp2)
  agilent_exp2 = apply(agilent_exp2,2,as.numeric)
  rownames(agilent_exp2) = gene

  #merge seperated agilent data
  agilent_exp2 = agilent_exp2[rownames(agilent_exp1),]
  table(rownames(agilent_exp1) == rownames(agilent_exp2))#T

  agilent_exp = cbind(agilent_exp1, agilent_exp2)

  #If NGS expression data is available, we prioritize its use.
  sample = setdiff(colnames(agliant_exp), colnames(rna_exp))
  m_agliant_exp = agliant_exp[,sample]
  
  #preprocessing
  rownames(gtex40) = gtex40$gene
  gtex40 = gtex40[-1,-1]
  pre_gtex40 = preprocess(gtex40)
  
  pre_m_agilent_exp = preprocess(m_agilent_exp)
  pre_rna_exp  = preprocess(rna_exp)
  
  #merge gene
  G = intersect(rownames(pre_m_agilent_exp), rownames(pre_rna_exp))
  G = intersect(rownames(pre_gtex40), G)
  
  pre_m_agilent_exp = pre_m_agilent_exp[G,]
  pre_rna_exp = pre_rna_exp[G,]
  pre_gtex40 = pre_gtex40[G,]
  
  #Final merge all files.
  test_gbm = cbind(pre_m_agilent_exp, pre_rna_exp, pre_gtex40)
}
