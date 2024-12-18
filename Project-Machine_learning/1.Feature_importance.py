import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.feature_selection import SelectFromModel
from joblib import Parallel, delayed
import sys


data = pd.read_csv('data/TCGA.integrative.model.input.by.rdgv.noSOD12.with.new_PPI.tsv',sep='\t') # Feature with PPI

# PPI related feature imputation by max value
data['avg_CPG'] =  data['avg_CPG'].fillna(np.max(data['avg_CPG']))
data['avg_DNA_repair'] =  data['avg_DNA_repair'].fillna(np.max(data['avg_DNA_repair']))
data['avg_PIK_mTOR'] =  data['avg_PIK_mTOR'].fillna(np.max(data['avg_PIK_mTOR']))
data['avg_Cell_Cycle'] =  data['avg_Cell_Cycle'].fillna(np.max(data['avg_Cell_Cycle']))
data['avg_ALL'] =  data['avg_ALL'].fillna(np.max(data['avg_ALL']))
data['distanceAll_path'] = data['distanceAll_path'].fillna(np.max(data['distanceAll_path'])) # just in case of distance All
data=data.drop(columns=['avg_CPG','avg_DNA_repair','avg_PIK_mTOR','avg_Cell_Cycle'])

# other Features imputation by zero
data=data.fillna(0) 

# import PPI features to choose one feature
# Feature_single = str(sys.argv[1]) 
Feature_single = 'distanceAll_path' # default => shortest distance without CPG


FF= Feature_single # save the file name with FF value

cpg_list= data[data['Is_positive']=='pos']['Gene'].values
other_list = data[data['Is_positive']=='neg']['Gene'].values

data = data[['Gene','case_con_zvalues', 'two_hit_zvalues', 'RDGV_by_nonRDGV', 'Tau_score',
       'oe_lof_upper', 'PTM_pvalue', 'PTV_frequency', 'amplification.freq',
       'deletion.freq', 'clin_zvalue', Feature_single,
       'Bidirectional_Deletion', 'Edition_efficiency', 'Guide_Abundance',
       'Insertion', 'Insertion_Deletion', 'PAM_Proximal_Deletion', 'Is_positive']]

repeat = 20000
splits = 5
INDEXXX = data.columns[1:-1]
np.random.seed(35)

NAME = ['RF_impu'] # model name to save -> randomforest imputation version used

# Parallelized loop using joblib 
# n
def process_iteration(i):
    selected_gene_other = np.random.choice(other_list, len(cpg_list), replace=False)
    fin_gene = np.concatenate((cpg_list, selected_gene_other))
    data3 = data[data.Gene.isin(fin_gene)]
    X = data3[data.columns[:-1]]
    y = data3[['Is_positive']]
    
    sss = StratifiedShuffleSplit(n_splits=splits, test_size=0.2, random_state=35)
    KKKK = []
    for j, (train, test) in enumerate(sss.split(X[X.columns[1:]], y['Is_positive'])):
        X_train_SS = X.iloc[train]
        y_train_SS = y.iloc[train]
        #X_test_SS = X.iloc[test]
        #y_test_SS = y.iloc[test]

        #before_fitting = time.time()
        sel = SelectFromModel(classifier)
        sel.fit(X_train_SS[X_train_SS.columns[1:]], y_train_SS['Is_positive'].values.ravel())


        importances = sel.estimator_.feature_importances_
        t_or_f = np.isin(sel.feature_names_in_, sel.get_feature_names_out()).astype(int)

        KKKK.append(importances)
        KKKK.append(t_or_f)
    return KKKK

# Create the RandomForestClassifier outside the loop
classifier = RandomForestClassifier(n_estimators=1000)
np.random.seed(35)
# Run the parallelized loop
results = Parallel(n_jobs=-1)(delayed(process_iteration)(i) for i in range(repeat)) # n_jobs =-1 might crash the cluster, reduce number based on the single sbatch job submission

# the way to save result file
reshaped_array = np.array(results).reshape(-1, np.array(results).shape[-1])
even_indices = reshaped_array[0::2]
odd_indices = reshaped_array[1::2]

even_df = pd.DataFrame(data=even_indices,columns=INDEXXX, ) # Even_df => the feature importance value by iteration
odd_df = pd.DataFrame(data=odd_indices,columns=INDEXXX, ) # Event of feature selection -> 1: selected 0: non-selected.

even_df.to_csv('./Result_feature_importance/%s.importance_PPI_%s.tsv.gz'%(NAME[0],FF ), sep='\t',compression='gzip')
odd_df.to_csv('./Result_feature_importance/%s.small_event_PPI_%s.tsv.gz'%(NAME[0],FF), sep='\t',compression='gzip')