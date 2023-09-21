#' @title SPde
#' @description
#' find Synergistic gene Pairs in Drug Efficacy.
#' Correlation analysis between CSGs and well-known RTK in GBM. 
#'
#' @param cox from output of CSGs.(gbm_agilent_ngs_gtex590_t_vs_n_cox.txt)
#' @param rtk; A data frame of well-known RTK genes.(RTK_list.csv)
#' @param exp; A data frame of TCGA-GBM merged expression. Extract only tumor data from output of MergeAgilentNGS using pheno data.
#'
#' @return A data.frame of filtered_corr.
#'
#' @export
SPde <- function(cox, rtk, exp){

  #Filter genes
  list_RTK_Cox = c(rtk$gene, cox)
  exp2 = exp[list_RTK_Cox,]
  
  # correlation analysis
  corr = data.frame(matrix(ncol=93, nrow=93))
  pval = data.frame(matrix(ncol=93, nrow=93))
  
  for (i in seq(91)){
    gene = list_RTK_Cox[i]
    for (j in seq(91)){
      cor_gene = list_RTK_Cox[j]
      correlation_result = cor.test(exp2[gene,], exp2[cor_gene,])
      cor = correlation_result$estimate
      p.val = correlation_result$p.value
      corr[i,j] = cor
      pval[i,j] = p.val
    }
  }
  
  # Filter results
  filtered_corr = corr
  filtered_corr[abs(corr) < 0.4] = NA
  filtered_corr[pval > 0.05] = NA
  
  return (filtered_corr)
}

