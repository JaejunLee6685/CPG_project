import numpy as np
import pandas as pd


common = [] 

for files in snakemake.input['vep_data_original'][0:]:
    k = pd.read_csv(files, sep='\t', low_memory=False, 
                usecols= ['ID','Gene','Feature','Feature_type','Consequence',
                          'Protein_position','cDNA_position' ,'CDS_position', 'Amino_acids','Codons',
                          'SYMBOL', 'SYMBOL_SOURCE',  'BIOTYPE', 'CANONICAL']).drop_duplicates()
    common.append(k)
com2 = pd.concat(common)
del(common)



rdgv = pd.read_csv(snakemake.input['rdgv'], sep='\t', compression='gzip', low_memory=False)


pposition = rdgv.merge(com2, left_on = ['ID','Gene','Feature','Feature_type','Consequence','SYMBOL', 'SYMBOL_SOURCE',  'BIOTYPE', 'CANONICAL'],
                       right_on= ['ID','Gene','Feature','Feature_type','Consequence','SYMBOL', 'SYMBOL_SOURCE',  'BIOTYPE', 'CANONICAL'],
                       
                       how='inner')
pposition['Protein_position'] = pposition['Protein_position'].astype(str)

pposition['PP_prot_pos'] = pposition['SYMBOL']+'_'+pposition['Protein_position']

pposition.to_csv(snakemake.output['pposition'],sep='\t', compression='gzip', index=False)