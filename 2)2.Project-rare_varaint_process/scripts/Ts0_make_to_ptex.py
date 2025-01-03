
import pandas as pd

data = pd.read_csv(snakemake.input['DATA'], 
                   compression='gzip', sep='\t', low_memory=False)
data2 = data[ (data['is_rdgv']==True) & (data['Consequence'].apply(lambda x: 'missense_variant' in x)) ]

data2.to_csv(snakemake.output['OUT'], 
             sep='\t', compression='gzip', index=False)