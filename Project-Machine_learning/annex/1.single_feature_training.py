# must include the way to calculate AUC, F1, Specificity, Sensitivity 
# based on AUC, permutation -> but the 

import numpy as np
import scipy.stats as sts
import pandas as pd

from sklearn.model_selection import  train_test_split
from sklearn.ensemble import RandomForestClassifier

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedShuffleSplit

from sklearn.metrics import roc_curve
from sklearn.metrics import auc

from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score

from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
import sys

data = pd.read_csv('data/TCGA.integrative.model.input.by.rdgv.noSOD12.with.new_PPI.tsv',
                   sep='\t',index_col=0)

# PPI related feature imputation by max value
data['avg_CPG'] =  data['avg_CPG'].fillna(np.max(data['avg_CPG']))
data['avg_DNA_repair'] =  data['avg_DNA_repair'].fillna(np.max(data['avg_DNA_repair']))
data['avg_PIK_mTOR'] =  data['avg_PIK_mTOR'].fillna(np.max(data['avg_PIK_mTOR']))
data['avg_Cell_Cycle'] =  data['avg_Cell_Cycle'].fillna(np.max(data['avg_Cell_Cycle']))
data['avg_ALL'] =  data['avg_ALL'].fillna(np.max(data['avg_ALL']))
data['distanceAll_path'] = data['distanceAll_path'].fillna(np.max(data['distanceAll_path']))

# other Features imputation by zero
data=data.fillna(0)


# import PPI features to choose one feature
#Feature_single = str(sys.argv[1])
Feature_single = 'distanceAll_path' # default => shortest distance without CPG


data = data[['case_con_zvalues', 'two_hit_zvalues', 'RDGV_by_nonRDGV', 'Tau_score',
       'oe_lof_upper', 'PTM_pvalue', 'PTV_frequency', 'amplification.freq',
       'deletion.freq', 'clin_zvalue', Feature_single,
       'Bidirectional_Deletion', 'Edition_efficiency', 'Guide_Abundance',
       'Insertion', 'Insertion_Deletion', 'PAM_Proximal_Deletion', 'Is_positive']] # Current version of selected features

# import features from sys argv when running sbatch
Feature_single = str(sys.argv[1]) # sys.argv[0] = python code itself / sys.argv[1] = parameter when submitting jobs
# example
# Feature_single = 'case_con_zvalues' 


# Feature_double = str(sys.argv[2])
# Feature_final = Feature_single+'_'+Feature_double

data = data[[Feature_single, 'Is_positive']] # select one feature (x) and y set (Y: positive or negative set)

# data = data[[Feature_single,Feature_double, 'Is_positive']]

name = 'RF_impu' # Machine learning model name to compare result

np.random.seed(35)
ML_pipe = RandomForestClassifier(n_estimators=1000,random_state=35) # random forest hyperparameter with default setting


np.random.seed(35)

first_permu = int(2e+4) # 20,000 repeation 
second_repeat = 5 # cross validation during the training model



x_col = [Feature_single, ]
#x_col = [Feature_single,Feature_double ]

data_y = data.pop('Is_positive') # pop y value
data_x = data # get x value (x value)


# normalization of dataset # to run machine learning model
scaler = StandardScaler()
data_for_training1 = pd.DataFrame(data = scaler.fit_transform(data_x ) ,index=data_x.index, 
            columns=data_x.columns)#.reset_index()

# new data with merged index
data3 = data_for_training1.merge(data_y, left_index=True, right_index=True)
data3['Is_positive'] = data3['Is_positive'].apply(lambda x: 1 if x =='pos' else 0)

# exclude TTN gene by manually
data3= data3[~data3.index.isin(['TTN'])]

cpg_list= data3[data3['Is_positive']==1].index.values
other_list = data3[data3['Is_positive']!=1].index.values


def process_iteration(i):
    DF_gene = []
    DF_score = []
    
    # selection of Machine learning set from Positive, negative sets
    selected_gene_other = np.random.choice(other_list, len(cpg_list), replace=False) 
    fin_gene = np.concatenate((cpg_list,selected_gene_other))
    data4= data3[data3.index.isin(fin_gene)].reset_index()
    
    # dataset division
    X= data4[np.concatenate((['Gene'],x_col))]
    y= data4[['Is_positive']]    
    sss= StratifiedShuffleSplit(n_splits=second_repeat, test_size=0.2)
    j=second_repeat-second_repeat
    
    FINAL_TABLE = []
    M2 = [] # AUC score save
    # repeat second_repeat times
    for train,test in sss.split(X[X.columns[1:]],y['Is_positive']) :

        X_train_SS = X.iloc[train]
        y_train_SS = y.iloc[train]
        X_test_SS = X.iloc[test]
        y_test_SS = y.iloc[test]
        
        # MODEL Fitting
        ML_pipe.fit(X_train_SS[X_train_SS.columns[1:]], y_train_SS['Is_positive'].values.ravel())
        
        # AUC
        ML_proba = ML_pipe.predict_proba(X_test_SS[X_test_SS.columns[1:]])
        malignant_probs = ML_proba[:,1]
        fpr, tpr, thresholds = roc_curve(y_test_SS, malignant_probs, pos_label=1)
        roc_auc = auc(fpr, tpr)
        auc_df = pd.DataFrame(data={'FPR':fpr, 'TPR':tpr, 'iters':'%s_%s'%(i,j),})
        M2.append(auc_df)        
        
        # PREDICTION and SCORE
        ML_pred = ML_pipe.predict(X_test_SS[X_test_SS.columns[1:]])
        precision = precision_score(y_test_SS, ML_pred)
        recall= recall_score(y_test_SS, ML_pred)
        f1score = f1_score(y_test_SS, ML_pred)
        accuracy =accuracy_score(y_test_SS, ML_pred)
        
        # table of likelihood
        daf = pd.DataFrame(data = {'Gene': data4[['Gene','Is_positive']].iloc[test].values[:,0],
                      'iters': np.repeat('%s_%s'%(i,j),len(test)),
                      'True_label': np.fromiter(data4[['Gene','Is_positive']].iloc[test].values[:,1],dtype=int),
                      'ML_label': np.fromiter(ML_pred,dtype=int) }
                          )
        # table of scores
        DDD = pd.DataFrame(data={'AUC':roc_auc, 'Precision':precision, "Recall":recall, 'F1_score':f1score,
                                'ACCURACY':accuracy}, index=['%s_%s'%(i,j)] )
        DF_gene.append(daf)
        DF_score.append(DDD)
        #SCORES.append([roc_auc, precision,recall,average_precision])
        j+=1
    FINAL_TABLE.append([DF_gene,DF_score,M2])
    #return DF_gene,DF_score
    return FINAL_TABLE   

result =  Parallel(n_jobs=-1)(delayed(process_iteration)(i) for i in range(first_permu)) # n_jobs = -1 use and ask all possible nodes from cluster / recommend not to use n=-1


re_result_shape =np.array(result,dtype='object').reshape(-1, np.array(result,dtype='object').shape[-1])

GENE = pd.concat([pd.concat(d) for d in re_result_shape[0::3]])
SCORE = pd.concat( [pd.concat(d) for d in re_result_shape[1::3]])
AUC_val= pd.concat( [pd.concat(d) for d in re_result_shape[2::3]])

GENE2 = GENE 

FIN_GENE = GENE2.groupby('Gene').sum(numeric_only=True)
FIN_GENE['EVENT'] = GENE2.groupby('Gene').count().values[:,1]
FIN_GENE['Likelyhood_CPG'] = FIN_GENE['ML_label'] / GENE2.groupby('Gene').count().values[:,1]


# save file vary on feature_single name from PPI

# Likelihood result after permutation, 
# CPG predicted as 1 / OTHER predicted as 0 
# original CPG - 1 / original OTH -1
FIN_GENE.to_csv('Result_machine_learning/%s_.ALL.permutation.%s.PPI.%s.tsv' %(name, (first_permu*second_repeat), Feature_single), sep='\t')

# double features
#FIN_GENE.to_csv('Result_machine_learning/%s_.ALL.permutation.%s.PPI.%s.tsv' %(name, (first_permu*second_repeat), Feature_final), sep='\t')

# Machine learning performance save 
SCORE.to_csv('Result_machine_learning/%s.ALL.score.%s.PPI.%s.tsv'%(name, (first_permu*second_repeat), Feature_single) , sep='\t')

# double features
#SCORE.to_csv('Result_machine_learning/%s.ALL.score.%s.PPI.%s.tsv'%(name, (first_permu*second_repeat), Feature_final) , sep='\t')


# TPR / FPR save after each iteration steps
AUC_val.to_csv('Result_machine_learning/%s.ALL.TPR_FPR.score.%s.PPI.%s.tsv'%(name, (first_permu*second_repeat), Feature_single) , sep='\t')

# double features
#AUC_val.to_csv('Result_machine_learning/%s.ALL.TPR_FPR.score.%s.PPI.%s.tsv'%(name, (first_permu*second_repeat), Feature_final) , sep='\t')

