import numpy as np
import pandas as pd


def cal_PTM_distance(PTM, variant): # by gene
    
    PTM_result = []
    
    for index, row in variant.iterrows():
        token = '-' != row['Protein_position']
        if token:
            value = '-' in row['Protein_position'] #multiple position
            if value:
                pos_min, pos_max = row['Protein_position'].split('-')
                if pos_min == '?':
                    pos_min=pos_max
                if pos_max == '?':
                    pos_max=pos_min
                
                for i, k  in PTM.iterrows():
                    ptm = k['position']
                    min_dif = np.abs(int(pos_min)-int(ptm))
                    max_dif = np.abs(int(pos_max)-int(ptm))


                    # if  (min_dif == 0) | (max_dif == 0):
                    #     PTM_result.append([row['ID'],row['SYMBOL'] ,'exact', ptm,min_dif, max_dif,row['Protein_position']])
                    #     #print(row['ID'],row['Gene.ensGene'] ,'exact', ptm,min_dif, max_dif,row['Protein_position'])
                        
                    # elif ( (min_dif >=1) & (min_dif <=2) ) | ( (max_dif >=1) & (max_dif <=2) ):
                    #     PTM_result.append([row['ID'],row['SYMBOL'] ,'proximal', ptm,min_dif, max_dif,row['Protein_position']])
                    #     #print(row['ID'],row['Gene.ensGene'],'proximal',ptm, min_dif, max_dif,row['Protein_position'])
                        
                    if ( (min_dif >=0) & (min_dif <=2) ) | ( (max_dif >=0) & (max_dif <=2) ):
                        PTM_result.append([row['ID'],row['SYMBOL'] ,'PTM_yes', ptm,min_dif, max_dif,row['Protein_position']])
                        #print(row['ID'],row['Gene.ensGene'],'distal', ptm,min_dif, max_dif,row['Protein_position'])
                        
                    elif (min_dif >= 3) & (max_dif >= 3):
                        continue
                        #print(row['ID'],row['Gene.ensGene'],'no_relation', ptm,min_dif, max_dif,row['Protein_position'])
            else:
                pos = row['Protein_position']
                for i, k in PTM.iterrows():
                    ptm = k['position']
                    diff = np.abs(int(ptm) - int(pos)) 
                    # if  diff == 0:
                    #     PTM_result.append([row['ID'],row['SYMBOL'],'exact', ptm,diff,diff,pos])
                    #     #print(row['ID'],row['Gene.ensGene'],'exact', ptm,diff,diff,pos)
                    # elif (diff >=1) & (diff <=2):
                    #     PTM_result.append([row['ID'],row['SYMBOL'],'proximal', ptm,diff,diff,pos])
                    #     #print(row['ID'],row['Gene.ensGene'],'proximal',ptm, diff,diff,pos)
                    if (diff >=0) & (diff <=2):
                        PTM_result.append([row['ID'],row['SYMBOL'],'PTM_yes', ptm,diff,diff,pos])
                        #print(row['ID'],row['Gene.ensGene'],'distal',ptm,diff,diff,pos)
                    else:
                        continue
                        #print(row['ID'],row['Gene.ensGene'],'no_relation',ptm,diff,pos)
        else:
            continue

    DF_TPM = pd.DataFrame(data=PTM_result, columns=['ID','SYMBOL','PTM_type','PTM_position',
                                                    'min_diff','max_diff','variant_distance'])
    
    dd = variant[['ID','SYMBOL','Protein_position']].merge(DF_TPM[['ID','SYMBOL','PTM_type','PTM_position']],
                                      left_on=['ID','SYMBOL'], right_on=['ID','SYMBOL'], how='left')
    dd['PTM_type'] = dd['PTM_type'].fillna('PTM_no')
    dd['PTM_position'] = dd['PTM_position'].fillna(0).astype(int)
    
    return(dd)

rdgv_fin = pd.read_csv(snakemake.input['pposition'], sep='\t',compression='gzip',low_memory=False)
dbptm = pd.read_csv(snakemake.input['PTM_DB'], sep='\t',low_memory=False)

big =[]
for gene in np.unique(rdgv_fin['SYMBOL']):
    if gene in dbptm['GENE'].values:
        PTM = dbptm[dbptm['GENE'] == gene]
        variant = rdgv_fin[rdgv_fin['SYMBOL'] == gene]
        
        gene1 = cal_PTM_distance(PTM, variant)
        big.append(gene1)
    else:
        continue
PTM_dbptm = pd.concat(big)
PTM_dbptm.to_csv(snakemake.output['PTM_VARIANT_OUT'], index=False, compression='gzip', sep='\t')