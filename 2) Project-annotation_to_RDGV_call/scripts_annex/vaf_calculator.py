import numpy as np
import pandas as pd

df = pd.read_csv(snakemake.input['tab'], sep='\t', low_memory=False)
df['POS'] = df['POS'].astype('int64')

def calculate_vaf(ad, rd):
    if ad == '.':
        return '.'
    aa = ad.split(',')
    if len(aa) == 2:
        return np.round(int(aa[1]) / (int(aa[0]) + int(aa[1])), 2)
    else:
        return np.round(int(ad) / (int(rd) + int(ad)), 2)
    

df['VAF'] = df.apply(lambda row: calculate_vaf(row['AD'], row['RD']), axis=1)
df= df.sort_values(by=['POS','SAMPLE'], ascending=True)

df.to_csv(snakemake.output['vaf'], sep = '\t', index = False, compression = 'gzip')

