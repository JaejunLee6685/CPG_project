import numpy as np
import pandas as pd

vep = pd.read_csv(snakemake.input['vep'], sep='\t', low_memory=False, compression='gzip')



'''ID      POS     FILTER  AC      gnomAD_FILTER   SAMPLE  GT      CANCER  ANCESTRY        CENTER  GENDER  Gene    
Feature Feature_type    Consequence     SYMBOL  SYMBOL_SOURCE   BIOTYPE CANONICAL       SIFT    PolyPhen      
EXON    INTRON  AF      AFR_AF  AMR_AF  EAS_AF  EUR_AF  SAS_AF  gnomADe_AF      gnomADe_AFR_AF  gnomADe_AMR_AF  
gnomADe_ASJ_AF  gnomADe_EAS_AF  gnomADe_FIN_AF  gnomADe_NFE_AF  
gnomADe_OTH_AF  gnomADe_SAS_AFSpliceAI_pred   REVEL   ClinVar ClinVar_CLNSIG  ClinVar_CLNREVSTAT      ClinVar_CLNDN

later when you merge 1KG and TCGA -> just take important value -> eventually to make VCF of total cohort -> recalculate to get cohort AF)
'''

AF_list = ['gnomADe_AF', 'gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_EAS_AF','gnomADe_NFE_AF','gnomADe_SAS_AF','gnomADe_OTH_AF',
           'TCGA_AFR_AF','TCGA_AMR_AF','TCGA_EAS_AF','TCGA_EUR_AF','TCGA_NAN_AF','TCGA_PAS_AF','TCGA_AF']
vep[AF_list] =vep[AF_list].replace('-',np.nan).apply(pd.to_numeric,errors='coerce')

# vep['gnomADe_AF'] =vep['gnomADe_AF'].replace('-', 0).astype(float)
# vep['gnomADe_AFR_AF'] =vep['gnomADe_AFR_AF'].replace('-', 0).astype(float)
# vep['gnomADe_AMR_AF'] =vep['gnomADe_AMR_AF'].replace('-', 0).astype(float)
# vep['gnomADe_EAS_AF'] =vep['gnomADe_EAS_AF'].replace('-', 0).astype(float)
# vep['gnomADe_NFE_AF'] =vep['gnomADe_NFE_AF'].replace('-', 0).astype(float)
# vep['gnomADe_OTH_AF'] =vep['gnomADe_OTH_AF'].replace('-', 0).astype(float)
# vep['gnomADe_SAS_AF'] =vep['gnomADe_SAS_AF'].replace('-', 0).astype(float)
# vep['gnomADe_FIN_AF'] =vep['gnomADe_FIN_AF'].replace('-', 0).astype(float)
# vep['gnomADe_ASJ_AF'] =vep['gnomADe_ASJ_AF'].replace('-', 0).astype(float)

common_05 = vep [
    
    (vep['gnomADe_AF'] >=0.05) &(vep['gnomADe_AFR_AF'] >=0.05) &(vep['gnomADe_AMR_AF'] >=0.05) &
    (vep['gnomADe_EAS_AF'] >=0.05) &(vep['gnomADe_NFE_AF'] >=0.05) &(vep['gnomADe_OTH_AF'] >=0.05) &  
    
    (vep['TCGA_AFR_AF'] >=0.05)&(vep['TCGA_AMR_AF'] >=0.05)&(vep['TCGA_EAS_AF'] >=0.05)&(vep['TCGA_EUR_AF'] >=0.05)&
    (vep['TCGA_NAN_AF'] >=0.05)&(vep['TCGA_PAS_AF'] >=0.05)&(vep['TCGA_AF'] >=0.05)
              ]


# only common but TCGA AF only used -> 
# when do you want to make AF fold change among Center -> [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 5.0, True]

#df2 = df.merge(common, left_on=['ID'], right_on=['ID'], how='inner')

common_01 = vep [
    
    (vep['gnomADe_AF'] >=0.01) &(vep['gnomADe_AFR_AF'] >=0.01) &(vep['gnomADe_AMR_AF'] >=0.01) &
    (vep['gnomADe_EAS_AF'] >=0.01) &(vep['gnomADe_NFE_AF'] >=0.01) &(vep['gnomADe_OTH_AF'] >=0.01) &  
    
    (vep['TCGA_AFR_AF'] >=0.01)&(vep['TCGA_AMR_AF'] >=0.01)&(vep['TCGA_EAS_AF'] >=0.01)&(vep['TCGA_EUR_AF'] >=0.01)&
    (vep['TCGA_NAN_AF'] >=0.01)&(vep['TCGA_PAS_AF'] >=0.01)&(vep['TCGA_AF'] >=0.01)
              ]

common_05 = common_05[common_05['Consequence'].apply(lambda x: 'missense_variant' in x)]
common_01 = common_01[common_01['Consequence'].apply(lambda x: 'missense_variant' in x)]


common_05.to_csv(snakemake.output['commons_05'], sep = '\t', index = False, compression = 'gzip')
common_01.to_csv(snakemake.output['commons_01'], sep = '\t', index = False, compression = 'gzip')


# # all temporary 
# df3 =df.merge(vep, left_on=['ID'], right_on=['ID'], how='inner' )
# df3.to_csv(snakemake.output['temporarys'], sep = '\t', index = False, compression = 'gzip')