import numpy as np
import pandas as pd


PTM_dbptm = pd.read_csv(snakemake.input['PTM_VARIANT_OUT'],  compression='gzip', sep='\t', low_memory=False)

rdgv_fin = pd.read_csv(snakemake.input['pposition'], sep='\t',compression='gzip',low_memory=False)

variant_all = pd.read_csv(snakemake.input['all_variant'], sep='\t', compression='gzip', low_memory=False)


# 1. /storage/scratch01/groups/gcc/LJJ_ML_project_2024/006.TCGA_VCF_FILTRATION/011.TCGA_rare_filter/006.get_dn_ds/out_merge/TCGA.gene.types.0.001.merge.tsv.gz

PTM_dbptm2 = PTM_dbptm.drop_duplicates(subset=['ID','SYMBOL','Protein_position','PTM_type'])

ptm_rdgv = PTM_dbptm2.merge(rdgv_fin,left_on=['ID','SYMBOL','Protein_position'],right_on=['ID','SYMBOL','Protein_position'], how='inner')

ptm_rdgv = ptm_rdgv.drop(columns=['SAMPLE','is_types'])

Merge_long = variant_all.merge(ptm_rdgv, left_on=['ID','SYMBOL'], right_on=['ID','SYMBOL'], how='left')

#groups = Merge_long.groupby(['SAMPLE','SYMBOL','is_types','PTM_type','ID']).count()



Merge_long2 = Merge_long.replace({"RV":"RDGV_no", "DGV":'RDGV_no', 'NV':'RDGV_no','RDGV':'RDGV_yes'})
# PTM_yes / PTM_no


Merge_long2_drop = Merge_long2.drop(columns=['SAMPLE'])
Merge_long2_drop =Merge_long2_drop.drop_duplicates()
M_fin = Merge_long2_drop[['ID','SYMBOL','is_types','PTM_type']]
M_fin.to_csv(snakemake.output['dropped_col'],sep='\t', compression='gzip',index=False)


Merge_long2_all = Merge_long2[['ID','SYMBOL','is_types','PTM_type']]
Merge_long2_all.to_csv(snakemake.output['nodrop_col'],sep='\t', compression='gzip',index=False)



# # 2. split them into RDGV, PTM

# rdgv_view = pd.pivot_table(Merge_long.groupby(['SAMPLE','SYMBOL','is_types']).count(),
#                            values='ID', columns = ['SYMBOL','is_types'], index=['SAMPLE'])

# rdgv_view_columns = rdgv_view.loc[:, rdgv_view.columns.get_level_values('is_types')=='RDGV']
# rdgv_view_columns.columns = [f"{gene}_RDGV" for gene in rdgv_view_columns.columns.get_level_values('SYMBOL')]
# rdgv_final = rdgv_view_columns.fillna(0).astype(int).applymap(lambda x: 1 if x>=1 else x )

# rdgv_final = rdgv_final[sorted(rdgv_final.columns)]

# ptm_view = pd.pivot_table(Merge_long.groupby(['SAMPLE','SYMBOL','PTM_type']).count(), 
#                           values='ID', columns = ['SYMBOL','PTM_type'], index=['SAMPLE'])

# ptm_view_columns = ptm_view.loc[:, ptm_view.columns.get_level_values('PTM_type')=='PTM_yes']
# ptm_view_columns.columns = [f"{gene}_PTM" for gene in ptm_view_columns.columns.get_level_values('SYMBOL')]
# ptm_final = ptm_view_columns.fillna(0).astype(int).applymap(lambda x: 1 if x>=1 else x )

# ptm_final = ptm_final[sorted(ptm_final.columns)]

# both_merge = pd.concat([rdgv_final,ptm_final], axis=1)

# #both_merge = both_merge[sorted(both_merge.columns)]

# # 3. to save the file 

# both_merge.to_csv(snakemake.output['RDGV_PTM'], sep='\t', compression='gzip')