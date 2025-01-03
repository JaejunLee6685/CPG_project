import pandas as pd
import numpy as np

longfile = pd.read_csv(snakemake.input['long'],sep='\t', compression='gzip', low_memory=False)
longfile44 =longfile
longfile =longfile.drop_duplicates(subset=['ID','FILTER','gnomAD_FILTER','SAMPLE'],keep=False) # different variant in the same gene

new_data = longfile.groupby(['ID','FILTER','CENTER','GT']).count()

new2 = pd.pivot_table(new_data, values='SAMPLE', index=['ID'], 
                      columns=['CENTER','GT'], 
            aggfunc=lambda x: int(x.values[0]) if len(x) > 0 else int(0)).fillna(0).astype(int).reset_index()

new2['BI_AL'] = new2[('BI')]['0/1'] + 2*new2[('BI')]['1/1'] 
new2['WUSM_AL'] = new2[('WUSM')]['0/1'] + 2*new2[('WUSM')]['1/1'] 
new2['BCM_AL'] = new2[('BCM')]['0/1'] + 2*new2[('BCM')]['1/1'] 
new2['WTSI_AL'] = new2[('WTSI')]['0/1'] + 2*new2[('WTSI')]['1/1'] 

new_data2 = new2[['ID', 'BI_AL', 'WUSM_AL', 'BCM_AL', 'WTSI_AL']]
new_data2['CENTER_AL'] = new_data2[['ID', 'BI_AL', 'WUSM_AL', 'BCM_AL', 'WTSI_AL']].sum(axis=1, numeric_only=True)

new_data2.loc[new_data2.index,'BCM_AF'] = (new_data2[['BCM_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['CENTER']).size()[0]) ).values
new_data2.loc[new_data2.index,'BI_AF'] = (new_data2[['BI_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['CENTER']).size()[1]) ).values
new_data2.loc[new_data2.index,'WTSI_AF'] = (new_data2[['WTSI_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['CENTER']).size()[2]) ).values
new_data2.loc[new_data2.index,'WUSM_AF'] = (new_data2[['WUSM_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['CENTER']).size()[3]) ).values
new_data2.loc[new_data2.index,'CENTER_AF'] = (new_data2[['CENTER_AL']]/ (2*longfile.drop_duplicates(['SAMPLE']).groupby(['CENTER']).size().sum()) ).values
new_data2.columns = new_data2.columns.get_level_values(0)

# ANCESTRY

wang_data = longfile.groupby(['ID','FILTER','ANCESTRY','GT']).count()

wanew2 = pd.pivot_table(wang_data, values='SAMPLE', index=['ID'], 
                      columns=['ANCESTRY','GT'], 
            aggfunc=lambda x: int(x.values[0]) if len(x) > 0 else int(0)).fillna(0).astype(int).reset_index()

wanew2['AFR_AL'] = wanew2 [('AFR')]['0/1'] + 2*wanew2 [('AFR')]['1/1'] 
wanew2['AMR_AL'] = wanew2[('AMR')]['0/1'] + 2*wanew2[('AMR')]['1/1'] 
wanew2['EAS_AL'] = wanew2[('EAS')]['0/1'] + 2*wanew2[('EAS')]['1/1'] 
wanew2['EUR_AL'] = wanew2[('EUR')]['0/1'] + 2*wanew2[('EUR')]['1/1'] 
wanew2['NAN_AL'] = wanew2[('NAN')]['0/1'] + 2*wanew2[('NAN')]['1/1'] 
wanew2['PAS_AL'] = wanew2[('PAS')]['0/1'] + 2*wanew2[('PAS')]['1/1'] 


wanew_data2 = wanew2[['ID', 'AFR_AL', 'AMR_AL', 'EAS_AL', 'EUR_AL','NAN_AL','PAS_AL']]

wanew_data2['TCGA_AL'] = wanew_data2[['ID', 'AFR_AL', 'AMR_AL', 'EAS_AL', 'EUR_AL','NAN_AL','PAS_AL']].sum(axis=1, numeric_only=True)

wanew_data2.loc[wanew_data2.index,'TCGA_AFR_AF'] = (wanew_data2[['AFR_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['ANCESTRY']).size()[0]) ).values
wanew_data2.loc[wanew_data2.index,'TCGA_AMR_AF'] = (wanew_data2[['AMR_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['ANCESTRY']).size()[1]) ).values
wanew_data2.loc[wanew_data2.index,'TCGA_EAS_AF'] = (wanew_data2[['EAS_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['ANCESTRY']).size()[2]) ).values
wanew_data2.loc[wanew_data2.index,'TCGA_EUR_AF'] = (wanew_data2[['EUR_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['ANCESTRY']).size()[3]) ).values
wanew_data2.loc[wanew_data2.index,'TCGA_NAN_AF'] = (wanew_data2[['NAN_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['ANCESTRY']).size()[4]) ).values
wanew_data2.loc[wanew_data2.index,'TCGA_PAS_AF'] = (wanew_data2[['PAS_AL']] / (2*longfile.drop_duplicates(['SAMPLE']).groupby(['ANCESTRY']).size()[5]) ).values

wanew_data2.loc[wanew_data2.index,'TCGA_AF'] = (wanew_data2[['TCGA_AL']]/ (2*longfile.drop_duplicates(['SAMPLE']).groupby(['ANCESTRY']).size().sum()) ).values

wanew_data2.columns = wanew_data2.columns.get_level_values(0)

king = new_data2.merge(wanew_data2, left_on = ['ID'], right_on =["ID"], how='inner')



again =longfile44

gg2 = again.merge(king, left_on=['ID'],right_on=['ID'],how='inner').sort_values(['ID','POS','SAMPLE'])

gg2.to_csv(snakemake.output['AF_out'], sep = '\t', index = False, compression = 'gzip',float_format='%.4f')