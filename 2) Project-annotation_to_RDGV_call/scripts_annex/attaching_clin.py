import numpy as np
import pandas as pd

df = pd.read_csv(snakemake.input['long'], sep='\t', low_memory=False, compression='gzip')

clin =pd.read_table(snakemake.input['clin'])

df2 = df.merge(clin[['SAMPLE','CANCER','ANCESTRY','CENTER','GENDER']], left_on =['SAMPLE'], right_on=['SAMPLE'], how='inner')
del(df)

df2=df2.sort_values(by=['POS','SAMPLE'], ascending=True)

df2.to_csv(snakemake.output['attached'], sep = '\t', index = False, compression = 'gzip')

