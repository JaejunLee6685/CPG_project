
import pandas as pd

data = pd.read_csv(snakemake.input['merged_all'],
                   compression='gzip', sep='\t', low_memory=False)
data2 = data[ (data['is_rdgv']==True) & (data['Consequence'].apply(lambda x: 'missense_variant' in x)) ]

data2.to_csv(snakemake.output['need2check']
             ,sep='\t', compression='gzip', index=False)