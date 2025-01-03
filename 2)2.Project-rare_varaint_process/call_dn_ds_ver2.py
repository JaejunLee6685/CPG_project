import numpy as np
import pandas as pd

col_to_read = ['ID','SAMPLE','CANCER','Gene','Feature','Feature_type','Consequence','SYMBOL', 'SYMBOL_SOURCE',  'BIOTYPE', 'CANONICAL','EXON','INTRON', ]

rare_damaging = pd.read_csv(snakemake.input['rare_damaging'],sep='\t', compression='gzip', low_memory=False, 
                            usecols = col_to_read)

non_rare_damgaing= pd.read_csv(snakemake.input['non_rare_damgaing'],sep='\t', compression='gzip', low_memory=False, 
                            usecols = col_to_read)

rare_non_damaging= pd.read_csv(snakemake.input['rare_non_damaging'],sep='\t', compression='gzip', low_memory=False, 
                            usecols = col_to_read)

non_rare_non_damaging= pd.read_csv(snakemake.input['non_rare_non_damaging'],sep='\t', compression='gzip', low_memory=False, 
                            usecols = col_to_read)


# make column to get types

rare_damaging['is_types'] = 'RDGV'
non_rare_damgaing['is_types'] = 'DGV'
rare_non_damaging['is_types'] = 'RV'
non_rare_non_damaging['is_types'] = 'NV'


merge = pd.concat([rare_damaging,non_rare_damgaing,rare_non_damaging,non_rare_non_damaging])

merge = merge.drop_duplicates(subset=['ID'])

merge.to_csv(snakemake.output['merge'],sep='\t', compression='gzip', index=False)

a_df = merge.groupby(['SYMBOL','is_types']).count()

b_df = pd.pivot_table(data=a_df, index=['SYMBOL'], columns=['is_types'], values=['ID'] ).fillna(0).astype(int)

b_df.to_csv(snakemake.output['count_gene_by_effect'], sep='\t', compression='gzip')