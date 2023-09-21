# SPde Main Functions
#
# Last revised date:
#   September 21, 2023
#
# Authors:
#   Juyoung Lee (wnddl111@catholic.ac.kr)
#   Hong lab team, catholic univ.

#' @title SPde
#' @description
#' Find Synergistic gene Pairs in Drug Efficacy 
#'
#' @param eff from output of PrepareCERES
#' @param annot from output of PrepareAnnotationFile
#' @param exp from output of MergeAgilentNGS
#' @param pheno; phenotype data(TCGA-GBM)
#' @param clinical; clincal data(TCGA-GBM)
#' @param cancer_type; The default is 'GBM', but you can replace it with another cancer type according to your research area.
#' @param gef; It means gene effect score. The default is -1, but you can replace it with another value.
#' @param percent; It represents what extent the distribution of GBM cell lines that meet the criteria will be. The default is 100%, but you can replace it with another value.
#' @param fc; It means fold change. The default is 1, but you can replace it with another value.
#'
#' @return A data.frame of res. This data frame represents the clinically significant genes in GBM that resulted from the Cox analysis.
#'
#' @export
SPde <- function(eff,annot,exp, pheno, clincal,cancer_type = 'GBM', gef = -1, percent = 1, fc = 1) {
  
	#1. Cell line level analysis
  ## Get GBM CERES data 
  sample <- annot[which(annot$Annotation == cancer_type), 'Depmap_ID']
  df <- eff[, sample]
  
  ## Count genes meeting the criteria
  cnt.table <- rowSums(df < gef)
  eGenes <- rownames(df[which(cnt.table >= percent * ncol(df)), ])
  
  message(paste0(cancer_type, 'cell line level analysis done..'))
  
  ## Extract genes which are important in survival of GBM cell lines
  eG <- intersect(rownames(exp), eGenes)
  exp <- exp[eG, ]
  
	#2. Patient level analysis
  ## Perform DEG analysis with extracted genes
  dres <- deg(exp, pheno, fc )
  new_genes <- dres[2][[1]]
  
  if (length(new_genes) <= 1) {
    message('no new gene..')
  } 
	else {
   c_exp <- as.data.frame(dres[1])
   colnames(c_exp) <- dres[3][[1]]
    
   message(paste0(cancer_type, ' DEG done..'))
   message(paste0(cancer_type, ' patient level analysis done..'))
    
   ## Read clinical data  
   mclinical <- clinical[, c('_PATIENT', 'OS', 'OS.time')]
   mclinical <- mclinical[!duplicated(mclinical), ]
    
   patient <- sapply(colnames(c_exp), function(x) substr(x, 1, nchar(x) - 3))
   t_c_exp <- as.data.frame(t(c_exp))
    
   ## Aggregate gene expression data
   m_final_vst <- aggregate(t_c_exp, list(patient), mean, na.rm = T)
   rownames(m_final_vst) <- m_final_vst$Group.1
   m_final_vst <- m_final_vst[, -1]
    
   analysis_gene <- intersect(new_genes, colnames(m_final_vst))
   g_m_final_vst <- m_final_vst[, analysis_gene]
   g_m_final_vst$'_PATIENT' <- rownames(g_m_final_vst)
    
   m_vst_clinical <- merge(g_m_final_vst, mclinical, by = '_PATIENT')
   rownames(m_vst_clinical) <- m_vst_clinical$'_PATIENT'
    
   ## Cox regression and survival analysis
   end = length(g_m_final_vst)
   covariates <- colnames(m_vst_clinical)[2:end]
  
   covariates <- lapply(covariates, function(x) gsub('-','.',x))
   colnames(m_vst_clinical)[2:end] = covariates
  
   univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(as.numeric(OS.time),OS)~',x)))
  
   univ_models <- lapply(univ_formulas,
                        function(x){coxph(x, data=m_vst_clinical)})
  
	   
   univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           p.value <- signif(x$wald['pvalue'],digits=2)
                           wald.test <- signif(x$wald['test'], digits = 2)
                           beta <- signif(x$coef[1], digits = 2)
                           HR <- signif(x$coef[2], digits = 2)
                           HR.confint.lower <- signif(x$conf.int[,'lower .95'],2)
                           HR.confint.upper <- signif(x$conf.int[,'upper .95'],2)
                           HR <- paste0(HR,' (',
                                        HR.confint.lower,'-',HR.confint.upper, ')')
                           res <- c(beta, HR, wald.test, p.value)
                           names(res) <-c('beta',' HR (95% CI for HR)', 'wald.test','p.value')
                           return(res)
                           
                         })
  
   res <- as.data.frame(t(as.data.frame(univ_results)))
   colnames(res) <- c('beta',' HR (95% CI for HR)', 'wald.test','p.value')
   rownames(res) <- covariates
  
	res$p.value = as.numeric(res$p.value)
  res_pval0.05 = res[ which(res$p.value <= 0.05 ),]
  storage.mode(res_pval0.05$beta) ='numeric'
  
  res_pval0.05$gene = rownames(res_pval0.05)
  res_pval0.05$col = sapply(res_pval0.05$beta, function(x) ifelse(x>0,'Positive','Negative'))
 }
  
  else{
    message(c,' no cox genes')
  }
return (res_pval0.05)
}
