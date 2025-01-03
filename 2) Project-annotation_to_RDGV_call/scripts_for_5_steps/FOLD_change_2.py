import numpy as np
import pandas as pd

common = [] 

for files in snakemake.input['commons'][0:]:
    k = pd.read_csv(files, sep='\t',compression='gzip', low_memory=False)
    common.append(k)
com2 = pd.concat(common)
del(common)



df = com2

df['AF_MAX'] = df[['BCM_KIT_AF','BI_KIT_AF','SeqCap2_AF','SeqCap3_AF','WTSI_KIT_AF','WUSM_KIT38_AF','WUSM_KIT18_AF']].max(axis=1)
df['AF_MIN'] = df[['BCM_KIT_AF','BI_KIT_AF','SeqCap2_AF','SeqCap3_AF','WTSI_KIT_AF','WUSM_KIT38_AF','WUSM_KIT18_AF']].min(axis=1)
df['Fold_change']  = np.around((df['AF_MAX']/df['AF_MIN']), 4)




# fold change based 

# for k,i in enumerate([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 5.0]):
#     data3 = df[
#         (df['Fold_change'] <= i) & (df['Fold_change'] !=np.inf)
#         ]
#     data3.to_csv('/storage/scratch01/groups/gcc/LJJ_ML_project_2024/006.TCGA_VCF_FILTRATION/008.get_exonic_splicing/out/s204_AF_FOLD/VEP.TCGA.COMMON.fold_%s_.tsv.gz'%k, 
#                  sep='\t',index = False, compression = 'gzip' )
#     print(len(data3.ID.unique()))



cuts = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 5.0]
if int(snakemake.params['folds2']) == 13:
    data2 = df
else:
    data2 = df[
        (df['Fold_change'] <= cuts[int(snakemake.params['folds2'])]) & (df['Fold_change'] !=np.inf)
        ]



data2.to_csv(snakemake.output['foldchange'],sep='\t', compression='gzip',index=False)