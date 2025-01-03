import numpy as np
import pandas as pd

data2 = pd.read_csv(snakemake.input['datas'],sep='\t', compression='gzip',low_memory=False)

data3 = pd.pivot_table(data2[['ID','gnomAD_FILTER', 'SAMPLE','GT','Consequence']],
                       values='GT', index=['ID','gnomAD_FILTER'], columns=['SAMPLE'],
                       aggfunc = lambda x: str(x.values[0]) if len(x) > 0 else '0/0')
data4 = data3.fillna('0/0')

col_sample = list(data2.SAMPLE.unique())

col_sample = list(np.sort(data2.SAMPLE.unique(), kind='quicksort') )

index1 = data4.reset_index()
del(data4)

index2 = index1
del(index1)

index2.loc[index2.index, '#CHROM'] = index2['ID'].str.split('_').str[0]
index2.loc[index2.index, 'POS'] = index2['ID'].str.split('_').str[1]
index2.loc[index2.index, 'REF'] = index2['ID'].str.split('_').str[2]
index2.loc[index2.index, 'ALT'] = index2['ID'].str.split('_').str[3]
index2.loc[index2.index, 'QUAL'] = '.'
index2.loc[index2.index, 'FILTER'] = 'PASS'
index2.loc[index2.index, 'INFO'] = '.'
index2.loc[index2.index, 'FORMAT'] = 'GT'

col_vcf = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO', 'FORMAT']
index3 = index2[col_vcf+col_sample]
index4 = index3.astype({'POS':'int'})

index4.sort_values(['#CHROM','POS']).to_csv(snakemake.output['vcf_tsv'], compression='gzip', sep='\t', index=0)

