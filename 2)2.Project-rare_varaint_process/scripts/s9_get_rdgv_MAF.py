import numpy as np
import pandas as pd


all_exon = pd.read_csv(snakemake.input['merged_all'],sep='\t', compression='gzip', low_memory=False)

# get only rdgv
all_exon_rdgv = all_exon[all_exon['is_rdgv']==True]
all_exon_rdgv = all_exon_rdgv.drop_duplicates()

# ptv determination
ptv = pd.read_csv(snakemake.input['merged_ptv'],sep='\t',compression='gzip',low_memory=False)
ptv =ptv.drop_duplicates()


ptv_pass = ptv[ptv['is_last_exon_intron']==False][['ID','SAMPLE','Gene']]
ptv_pass = ptv[ptv['is_last_exon_intron']==False][['ID','SAMPLE','Gene']].drop_duplicates()

filter_by_ptv = all_exon_rdgv[all_exon_rdgv[['ID', 'SAMPLE','Gene']].apply(tuple, axis=1).isin(ptv_pass[['ID', 'SAMPLE','Gene']].apply(tuple,axis=1))]


# pext filtration
ptex = pd.read_csv(snakemake.input['pext_fin'], sep='\t', compression='gzip',low_memory=False)
ptex_pass = ptex[ptex['is_ptex']==True]

filter_by_ptex = all_exon_rdgv[all_exon_rdgv[['ID', 'SAMPLE','Gene']].apply(tuple, axis=1).isin(ptex_pass[['ID', 'SAMPLE','Gene']].apply(tuple,axis=1))]

# other filtration

filter_by_ohter= all_exon_rdgv[(all_exon_rdgv['is_clin']==True) & 
               ~(all_exon_rdgv['Consequence'].apply(lambda x: 'missense_variant' in x)) &
              (all_exon_rdgv['is_ptv']!=True) ] # other


FINAL = pd.concat([filter_by_ptv,filter_by_ptex,filter_by_ohter ])
FINAL = FINAL.drop_duplicates()

FINAL.to_csv(snakemake.output['RDGV'], sep='\t', compression='gzip', index=False)