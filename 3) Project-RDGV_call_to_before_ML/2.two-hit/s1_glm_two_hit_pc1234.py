import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import dmatrices
import scipy.stats as stats


# import clinical data
clin =pd.read_csv(snakemake.input['clin'],sep='\t', 
                  usecols=['SAMPLE_ID_ONLY','CANCER','ANCESTRY','CENTER','GENDER',
                          'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10',])

# TCGA only clinical data -> SAMPLE_ID_ONLY -to merge germline and somatic

# import germline data and change sampleID
germline = pd.read_csv(snakemake.input['germline'], sep='\t', compression='gzip', 
                       low_memory = False, )
germline['SAMPLE_raw'] = germline['SAMPLE'].apply(lambda x: x[0:12])
germline.index = germline['SAMPLE_raw']
germline = germline.drop(columns=['SAMPLE_raw'])
germline = germline.rename(columns ={'SAMPLE':"SAMPLE_ID"})

# import somatic loh 
loh = pd.read_csv(snakemake.input['somatic'],sep='\t', compression='gzip',         # './TCGA.LOH.ALL.tsv.gz'
                          low_memory=False)
loh['SAMPLE_raw'] = loh['SAMPLE_ID'].apply(lambda x: x[0:12])
loh.index = loh['SAMPLE_raw']
loh = loh.drop(columns=['SAMPLE_raw'])

# matching sample
germ_somatic_matching = germline.index.intersection(loh.index)

# make new genes  - germline
germline_sample_gene = germline[loh.columns.intersection(germline.columns)].loc[germ_somatic_matching]


# make new genes - somatic
loh_sample_gene = loh[loh.columns.intersection(germline.columns)].loc[germ_somatic_matching]

loh_sample_gene =pd.DataFrame(data=loh_sample_gene.values, index=loh_sample_gene.index, 
                              columns=[i+'_LOH' for i in list(loh_sample_gene.columns)] ) # make LOH comments -for N in glm
loh_sample_gene =loh_sample_gene.rename(columns={'SAMPLE_ID_LOH':'SAMPLE_ID'}) 

# make input 
merge_loh_germline = germline_sample_gene.merge(loh_sample_gene, left_index=True, right_index=True)
merge_loh_germline = merge_loh_germline.drop(columns=['SAMPLE_ID_y'])
merge_loh_germline = merge_loh_germline.rename(columns={'SAMPLE_ID_x':'SAMPLE_ID'}) 
merge_loh_germline2 = merge_loh_germline.drop(columns=['SAMPLE_ID'])
del(merge_loh_germline) # required?

merge_loh_germline2 =merge_loh_germline2.astype(int)


# make INPUT with pc
INPUT_ver1 = merge_loh_germline2.merge(clin[['SAMPLE_ID_ONLY','CANCER',
                            'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10' ]], 
                          left_index=True, right_on=['SAMPLE_ID_ONLY'])
INPUT_ver1.index = INPUT_ver1.SAMPLE_ID_ONLY
INPUT_ver1 = INPUT_ver1.drop(columns=['SAMPLE_ID_ONLY'])
INPUT_ver1.columns = INPUT_ver1.columns.str.replace('-','__') # change multipe gene-gene

# change the order of columns -11 -> Cancer + PC1~10 
concat_gene_pc_cancer = np.concatenate((
    INPUT_ver1.columns[-11:] , 
    INPUT_ver1.columns[:-11]
              ))

# input version two with re-order
INPUT_ver2 = INPUT_ver1[concat_gene_pc_cancer]

# input columns order 
# GENE
first_col = merge_loh_germline2.columns[0:int(len(merge_loh_germline2.columns)/2)]
first_col = first_col.str.replace('-','__')

# GENE_LOH
second_col = merge_loh_germline2.columns[int(len(merge_loh_germline2.columns)/2):]
second_col = second_col.str.replace('-','__') 

PCs = 'PC1'+ '+' + 'PC2' + '+' + 'PC3' + '+' + 'PC4'

# model construct 
big=[]
k=0
for i,j in zip(first_col, second_col): # j = GENE_LOH / i = GENE
    #print(i,j)
    ty, tX = dmatrices('%s ~ %s + %s' %(j, i, PCs), data=INPUT_ver2, return_type='dataframe')
    tmm3 = sm.GLM(ty,tX, family=sm.families.Binomial())
    result_mm3 =tmm3.fit()
    
    # add data
    coeffcient = result_mm3.params[1]
    std_err = result_mm3.bse[1]
    zvalues = result_mm3.tvalues[1]
    pvalues = result_mm3.pvalues[1]
    left= result_mm3.conf_int().iloc[1,:].values[0]
    right= result_mm3.conf_int().iloc[1,:].values[1]
    deviance =  result_mm3.deviance
    aic= result_mm3.aic
    bic=result_mm3.bic
    llf = result_mm3.llf
    p_chi2 = result_mm3.pearson_chi2
    p_r2 = result_mm3.pseudo_rsquared()
    
    big.append([coeffcient,std_err,zvalues,pvalues,left,right,deviance,aic,bic,llf,p_chi2,p_r2])
#     k+=1
#     if k >10:
#         break

bb = pd.DataFrame(data=big, index=merge_loh_germline2.columns[0:int(len(merge_loh_germline2.columns)/2)], 
                 columns = ['coeffcient', 'std_err', 'zvalues', 'pvalues', '[0.025]', '[0.0975]',
     'deviance', 'aic','bic', 'log-likelihood', 'pearsons chi2','pseudo_r2'])

bb.to_csv(snakemake.output['two_hit_glm'], sep='\t')