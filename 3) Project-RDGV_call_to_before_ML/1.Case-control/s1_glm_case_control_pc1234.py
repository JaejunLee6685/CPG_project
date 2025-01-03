import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import dmatrices
import scipy.stats as stats
import warnings
warnings.filterwarnings('ignore')

# case-control preparation

clin = pd.read_csv(snakemake.input['clin'], sep='\t', 
                   usecols=['SAMPLE','CANCER','ANCESTRY','CENTER','GENDER',
                          'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10',])

clin['CASE_CONTROL'] = clin['CANCER'].apply(lambda x: 0 if x=='CONTROL' else 1)
clin=clin[['SAMPLE','CASE_CONTROL','CANCER','ANCESTRY','CENTER','GENDER',
                          'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']]
#print(clin.shape)
clin=clin[clin['SAMPLE']!='TCGA-29-2429-10A-01D-1526-09'] # outlier from TCGA cohort
#print(clin.shape)

#######################

# import case -> depends on the MAF -> mutation file required
case = pd.read_csv(snakemake.input['cases'], sep='\t', compression='gzip',
                  low_memory=False, index_col=0)

# get left join to get all tcga samples 
clin_case = clin[clin['CASE_CONTROL']==1]
clin_case_mut = clin_case.merge(case, left_on=['SAMPLE'], right_index=True, how='left')
clin_case_mut=clin_case_mut.fillna(0)

#######################

# import control -> depends on the MAF
control = pd.read_csv(snakemake.input['control'], sep='\t', compression='gzip',
                     low_memory=False, index_col=0)

# get left join to get all 1kg samples
clin_control = clin[clin['CASE_CONTROL']==0]
clin_control_mut = clin_control.merge(control, left_on=['SAMPLE'], right_index=True, how='left')
clin_control_mut=clin_control_mut.fillna(0)

# match same gene from case-control
not_in_case_gene = clin_case_mut.columns.union(clin_control_mut.columns).difference(clin_case_mut.columns)
not_in_control_gene = clin_case_mut.columns.union(clin_control_mut.columns).difference(clin_control_mut.columns)

# um mapped genes -> 0 value
clin_case_mut[not_in_case_gene] = 0
clin_control_mut[not_in_control_gene] = 0

# make input for glm

INPUT_ver1 = pd.concat([clin_control_mut,clin_case_mut])

INPUT_ver1 =INPUT_ver1.sort_values(['CASE_CONTROL','SAMPLE'],ascending=[False,True] )
INPUT_ver1.columns = INPUT_ver1.columns.str.replace('-','__') # to change geneA-GeneB -> geneA__gene_B (if not regression model recognize geneA-GeneB as the arithmetic equation)

# INPUT_ver1.columns[17:]
# consider only genes (17 -> A1BG to ZZZ3 ) / the order will not change until I include more PC values


##############################
# model construction
# for loop based on genes -> takes a few moments

PCs = 'PC1'+ '+' + 'PC2' + '+' + 'PC3' + '+' + 'PC4'
# pcs can be changable - by comparing models' pvalue or others 

big=[]
k=0
for i in INPUT_ver1.columns[16:] :
    #print(i,j)
    ty, tX = dmatrices('CASE_CONTROL ~ %s + %s' %(i,PCs), data=INPUT_ver1, return_type='dataframe')
    tmm3 = sm.GLM(ty,tX, family=sm.families.Binomial())
    result_mm3 =tmm3.fit()
    
    # add data #by gene
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

# to generate dataframe
bb = pd.DataFrame(data=big, index=INPUT_ver1.columns[16:].str.replace('__','-'), 
                 columns = ['coeffcient', 'std_err', 'zvalues', 'pvalues', '[0.025]', '[0.0975]',
     'deviance', 'aic','bic', 'log-likelihood', 'pearsons chi2','pseudo_r2'])

bb.to_csv(snakemake.output['case_control_glm'], sep='\t')