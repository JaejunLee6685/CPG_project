import numpy as np
import pandas as pd
import scipy.stats as sts


def perform_fisher_test(gene, data_table):
    each_gene = data_table[data_table['SYMBOL']==gene]
    
    OO = each_gene[(each_gene['is_types']=='RDGV_yes') & (each_gene['PTM_type']=='PTM_yes')].shape[0]
    OX = each_gene[(each_gene['is_types']=='RDGV_yes') & (each_gene['PTM_type']!='PTM_yes')].shape[0]
    XO = each_gene[(each_gene['is_types']!='RDGV_yes') & (each_gene['PTM_type']=='PTM_yes')].shape[0]
    XX = each_gene[(each_gene['is_types']!='RDGV_yes') & (each_gene['PTM_type']!='PTM_yes')].shape[0]
    
    # exact distance
    exact_fisher = sts.fisher_exact([[OO,OX],[XO,XX]], alternative='greater')
    
    return ([gene, exact_fisher[0], exact_fisher[1]])


data = pd.read_csv(snakemake.input['PTM_distance'], sep='\t', compression='gzip', low_memory=False)

gene_list = np.unique(data['SYMBOL'])

big = []

fisher_test_all = list(map(lambda gene: perform_fisher_test(gene,data), gene_list))
df_of_fisher= pd.DataFrame(data=fisher_test_all,columns=['Gene','odd_ratio','pvalue']).set_index('Gene')
big.append(df_of_fisher)


big_df = pd.concat(big, axis=1)

big_df.to_csv(snakemake.output['Fisher'], sep='\t', )