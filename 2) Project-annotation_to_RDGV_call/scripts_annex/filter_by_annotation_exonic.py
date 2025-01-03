import numpy as np
import pandas as pd

df = pd.read_csv(snakemake.input['wide_long'], sep='\t', low_memory=False, compression='gzip')

exonic =pd.read_csv(snakemake.input['filter_ed'],sep='\t', low_memory=False)

df2 = df[df['ID'].isin(exonic['ID'])]
df2 =df2.sort_values(['POS','SAMPLE'], ascending=True)

df2.to_csv(snakemake.output['long_filter'], sep = '\t', index = False, compression = 'gzip')

