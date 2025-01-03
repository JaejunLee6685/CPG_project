import numpy as np
import pandas as pd

#df = pd.read_csv(snakemake.input['wide_long'], sep='\t', low_memory=False, compression='gzip')

df =pd.read_csv(snakemake.input['filter_ed'],sep='\t', low_memory=False, compression='gzip')

#df2 = df[df['ID'].isin(exonic['ID'])]
df2 =df.sort_values(['POS','SAMPLE'], ascending=True)

df3 = df2


# df2 = df2[df.duplicated(subset=['ID','FILTER','gnomAD_FILTER'], keep=False)]
# DD = df2[~df.duplicated(subset=['ID','FILTER','gnomAD_FILTER'], keep=False)]


# data2 = df2.drop_duplicates(['ID','FILTER','gnomAD_FILTER'])
# data3 =data2[data2.duplicated(subset=['ID'],keep=False)].sort_values(by=['ID','gnomAD_FILTER'], ascending = [True, False])


gnomad_filter_dict = {}

# First pass to collect the gnomAD_FILTER values
for index, row in df3.iterrows():
    if row['gnomAD_FILTER'] != '.':
        gnomad_filter_dict[row['ID']] = row['gnomAD_FILTER']

# Second pass to update the gnomAD_FILTER values
for index, row in df3.iterrows():
    if row['gnomAD_FILTER'] == '.' and row['ID'] in gnomad_filter_dict:
        df3.at[index, 'gnomAD_FILTER'] = gnomad_filter_dict[row['ID']]

# data4 = data3

# data5 = pd.concat([DD,data4])
# data5= data5.sort_values(['POS','SAMPLE'], ascending=True)

#df2.to_csv(snakemake.output['long_filter'], sep = '\t', index = False, compression = 'gzip')

df3.to_csv(snakemake.output['long_filter'], sep = '\t', index = False, compression = 'gzip')

#data5.to_csv(snakemake.output['long_filter'], sep = '\t', index = False, compression = 'gzip')