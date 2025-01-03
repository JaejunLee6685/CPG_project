import numpy as np
import pandas as pd

def get_second_value(cell):
    parts = cell.split(',')
    return int(parts[1]) if len(parts) > 1 else int(cell)

df = pd.read_csv(snakemake.input['wide_long'], sep='\t', low_memory=False, compression='gzip')

# drop no data in DP, VAF, GQ
df = df[(df['DP']!='.') & (df['VAF']!='.') &(df['GQ']!='.') ]
df['DP'] = df['DP'].astype(int)
df['VAF'] = df['VAF'].astype(float)
df['GQ'] = df['GQ'].astype(int)

# split with GT 01 11
# GT 0/1 -> GQ 20 ,DP 10, 0.2 VAF 0.8
df2 = df[ ( (df['gnomAD_FILTER'] =='PASS') | (df['gnomAD_FILTER'] =='.') ) & (df['GQ'] >=20) & (df['DP'] >=10) & 
         
        (
        ( (df['GT'] =="0/1" ) & (df['VAF'] >= 0.2) & (df['VAF'] <= 0.8) ) 
        )
         
         
         ]

# GT 1/1 -> GQ 20 ,DP 10, 0.8 VAF
df22 = df[ ( (df['gnomAD_FILTER'] =='PASS') | (df['gnomAD_FILTER'] =='.') ) & (df['GQ'] >=20) & (df['DP'] >=10) & 
         
        (
        ( (df['GT'] =="1/1" ) & (df['VAF'] > 0.8) )
        )
         
         
         ]
#DF1 = pd.concat([df2,df22]).sort_values(['POS','SAMPLE']) # AD not consider


df3 =df2

df3['ALT'] = df3['AD'].apply(get_second_value)
df4 =df3[df3['ALT']>=5]
df4 =df4.drop(columns=["ALT"])


DF2 = pd.concat([df4,df22]).sort_values(['POS','SAMPLE']) # AD >=5

data2 = DF2.drop_duplicates(['ID','FILTER','gnomAD_FILTER'])

dupee = data2[data2.duplicated(['ID'],keep=False)].sort_values('ID')

DF3 = DF2[~DF2.ID.isin(dupee.ID.unique())]

#DF1.to_csv(snakemake.output['DF1'], sep = '\t', index = False, compression = 'gzip')
DF3.to_csv(snakemake.output['DF3'], sep = '\t', index = False, compression = 'gzip')