# SPde

```SPde``` is an source code to discover potential synergistic pair genes that could modulate drug resistance mechanisms in GBM patients undergoing RTK targeted therapy. We hypothesized that such regulators could exert a pivotal influence on the overall survival of GBM patients. To achieve this, we conducted a two-stage analysis, first with cell line level analysis then patient level analysis. 


![image](https://github.com/Honglab-Research/SPde/assets/79962288/fccac2a3-3ec6-4b66-bfe7-6c50c64e2876)

## 01. Dependencies
- R (>= 4.0.0)
- forcats (1.0.0)
- survival (3.4-0)
- ggplot2 (3.4.2)
- dplyr (1.1.0)
- survminer (0.4.9)
- limma (3.54.2)
- survival (3.4-0)
- readxl (1.4.2)
  
## 02. Usage

All input files are in ```data``` folder. Please refer to the example scripts in the ```examples``` folder.

### A. Load data
Set the working directory to where you downloaded these input files and load the necessary data files.
```
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

```

### B. Merge TCGA-GBM RNAseq data
Merged separated NGS and Agilent TCGA-GBM data, and added GTEx data for DEG analysis.
```
library(MutaliskR)
library(BSgenome.Hsapiens.UCSC.hg38)

annovar.file <- system.file("extdata/03_Samples", "HCC1187--HCC-1187BL.snv.indel.final.v6.annotated.hg38_multianno.txt", package = "MutaliskR")
input <- PrepareAnnovarFile(annovar.file = annovar.file)
results <- IdentifyDbsSignatures(
  input = input,
  bsg = BSgenome.Hsapiens.UCSC.hg38,
  analyze.variants = TRUE,
  analyze.variants.column.gene = "Gene.refGene",
  analyze.variants.column.group = "Func.refGene",
  sample.id = "HCC1187",
  n.cores = 2,
  save.dir = "Output directory path"
)
```

### C. Preprocessing
Excluded genes with zero expression in all samples and performed normalization.
The source code for ```PrepareCERES``` and ```PrepareAnnotationFile``` can be found in the R folder.
```
eff <- PrepareCERES(eff)
annot <- PrepareAnnotationFile(annot)
```

### D. Run
You can obtain synergistic genes with well-known RTK in GBM using the following code.
Make sure to manually compare the prognosis difference between synergistic pair genes when they act together and when they act individually.
The source code for ```CSGs``` and ```SPde``` can be found in R folder.
```
#Run CSGs
csg <- CSGs(eff,annot,exp, pheno, clincal,cancer_type = 'GBM', gef = -1, percent = 1, fc = 1) 

#Run SPde
SP <- SPde(csg$gene, rtk, exp)
```


## 03. Contact
For any questions or issues, please feel free to contact us or post them in the Issues tab of this repository.
