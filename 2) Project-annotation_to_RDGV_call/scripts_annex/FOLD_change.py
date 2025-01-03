import numpy as np
import pandas as pd

common = [] 

for files in snakemake.input['commons'][0:]:
    k = pd.read_csv(files, sep='\t',compression='gzip', low_memory=False)
    common.append(k)
com2 = pd.concat(common)
del(common)



df = com2

df['AF_MAX'] = df[['BI_AF','WUSM_AF','BCM_AF','WTSI_AF']].max(axis=1)
df['AF_MIN'] = df[['BI_AF','WUSM_AF','BCM_AF','WTSI_AF']].min(axis=1)
df['Fold_change']  = np.around((df['AF_MAX']/df['AF_MIN']), 4)


# fold change based 

for k,i in enumerate([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 5.0,10]):
    if i == 10:
         data3.to_csv('out5/s004_AF_FOLDchange/VEP.TCGA.COMMON.fold_%s_.tsv.gz'%k, 
                 sep='\t',index = False, compression = 'gzip' ) # original foldchange -> raw VCF

    data3 = df[
        (df['Fold_change'] <= i) & (df['Fold_change'] !=np.inf)
        ]
    #/VEP.TCGA.COMMON.fold_{folls}_.tsv.gz
    data3.to_csv('out5/s004_AF_FOLDchange/VEP.TCGA.COMMON.fold_%s_.tsv.gz'%k, 
                 sep='\t',index = False, compression = 'gzip' )
    print(len(data3.ID.unique()))