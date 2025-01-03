import numpy as np
import pandas as pd

df = pd.read_csv(snakemake.input['DF3'], sep='\t', low_memory=False, compression='gzip')

vep =pd.read_table(snakemake.input['vep'],sep='\t',low_memory=False, 
                   usecols = ['ID','Gene','Feature','Feature_type','Consequence','SYMBOL',
                              'SYMBOL_SOURCE','BIOTYPE','CANONICAL','EXON', 'INTRON', 
                              'gnomADe_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_EAS_AF','gnomADe_NFE_AF','gnomADe_OTH_AF','gnomADe_SAS_AF','gnomADe_FIN_AF','gnomADe_ASJ_AF',
                              'SIFT','PolyPhen','SpliceAI_pred',  'REVEL',
                              'ClinVar', 'ClinVar_CLNSIG',  'ClinVar_CLNREVSTAT',   'ClinVar_CLNDN'
])

vep['gnomADe_AF'] =vep['gnomADe_AF'].replace('-', 0).astype(float)
vep['gnomADe_AFR_AF'] =vep['gnomADe_AFR_AF'].replace('-', 0).astype(float)
vep['gnomADe_AMR_AF'] =vep['gnomADe_AMR_AF'].replace('-', 0).astype(float)
vep['gnomADe_EAS_AF'] =vep['gnomADe_EAS_AF'].replace('-', 0).astype(float)
vep['gnomADe_NFE_AF'] =vep['gnomADe_NFE_AF'].replace('-', 0).astype(float)
vep['gnomADe_OTH_AF'] =vep['gnomADe_OTH_AF'].replace('-', 0).astype(float)
vep['gnomADe_SAS_AF'] =vep['gnomADe_SAS_AF'].replace('-', 0).astype(float)
vep['gnomADe_FIN_AF'] =vep['gnomADe_FIN_AF'].replace('-', 0).astype(float)
vep['gnomADe_ASJ_AF'] =vep['gnomADe_ASJ_AF'].replace('-', 0).astype(float)

common = vep [
    
    (vep['gnomADe_AF'] >=0.05) &(vep['gnomADe_AFR_AF'] >=0.05) &(vep['gnomADe_AMR_AF'] >=0.05) &
    (vep['gnomADe_EAS_AF'] >=0.05) &(vep['gnomADe_NFE_AF'] >=0.05) &(vep['gnomADe_OTH_AF'] >=0.05) 
              ]


# only common but TCGA AF only used -> 
# when do you want to make AF fold change among Center -> [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 5.0, True]

df2 = df.merge(common, left_on=['ID'], right_on=['ID'], how='inner')

df2 = df2[
    
    (df2['TCGA_AFR_AF'] >=0.05)&(df2['TCGA_AMR_AF'] >=0.05)&(df2['TCGA_EAS_AF'] >=0.05)&(df2['TCGA_EUR_AF'] >=0.05)&
    (df2['TCGA_NAN_AF'] >=0.05)&(df2['TCGA_PAS_AF'] >=0.05)&(df2['TCGA_AF'] >=0.05)

          ]


df25 = df2[df2['Consequence'].str.contains('miss',regex=False)]

df25.to_csv(snakemake.output['commons'], sep = '\t', index = False, compression = 'gzip')

# all temporary 
df3 =df.merge(vep, left_on=['ID'], right_on=['ID'], how='inner' )
df3.to_csv(snakemake.output['temporarys'], sep = '\t', index = False, compression = 'gzip')