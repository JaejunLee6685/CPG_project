import numpy as np
import pandas as pd

common = [] 

for files in snakemake.input['data'][0:]:
    k = pd.read_csv(files, sep='\t',compression='gzip', low_memory=False)
    common.append(k)
com2 = pd.concat(common)
del(common)


com2.sort_values(by=['ID','POS','SAMPLE']).to_csv(snakemake.output['merged'], compression='gzip', sep='\t', index=0)