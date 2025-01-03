from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedShuffleSplit, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, precision_score, recall_score, f1_score, accuracy_score
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

data = pd.read_csv('./TCGA.integrative.model.input.by.rdgv.noSOD12.with.new_PPI.tsv', sep='\t', index_col=0)

data['avg_CPG'] = data['avg_CPG'].fillna(np.max(data['avg_CPG']))
data['avg_DNA_repair'] = data['avg_DNA_repair'].fillna(np.max(data['avg_DNA_repair']))
data['avg_PIK_mTOR'] = data['avg_PIK_mTOR'].fillna(np.max(data['avg_PIK_mTOR']))
data['avg_Cell_Cycle'] = data['avg_Cell_Cycle'].fillna(np.max(data['avg_Cell_Cycle']))
data['avg_ALL'] = data['avg_ALL'].fillna(np.max(data['avg_ALL']))
data['distanceAll_path'] = data['distanceAll_path'].fillna(np.max(data['distanceAll_path']))
data = data.fillna(0)

data = data[['case_con_zvalues', 'two_hit_zvalues', 'RDGV_by_nonRDGV', 'Tau_score',
       'oe_lof_upper', 'PTM_pvalue', 'PTV_frequency', 'amplification.freq',
       'deletion.freq', 'clin_zvalue', 'distanceAll_path',
'Bidirectional_Deletion', 'Edition_efficiency', 'Guide_Abundance',
       'Insertion', 'Insertion_Deletion', 'PAM_Proximal_Deletion', 'Is_positive']]

np.random.seed(35)

name = 'XGBoost'
first_permu = int(2e+4)
second_repeat = 5

x_col = ['case_con_zvalues', 'two_hit_zvalues', 'RDGV_by_nonRDGV', 'Tau_score',
       'oe_lof_upper', 'PTM_pvalue', 'PTV_frequency', 'amplification.freq',
       'deletion.freq', 'clin_zvalue', 'distanceAll_path',
'Bidirectional_Deletion', 'Edition_efficiency', 'Guide_Abundance',
       'Insertion', 'Insertion_Deletion', 'PAM_Proximal_Deletion', ]

data_y = data.pop('Is_positive')
data_x = data

scaler = StandardScaler()
data_for_training1 = pd.DataFrame(data=scaler.fit_transform(data_x), 
                                  index=data_x.index, columns=data_x.columns)
data3 = data_for_training1.merge(data_y, left_index=True, right_index=True)
data3['Is_positive'] = data3['Is_positive'].apply(lambda x: 1 if x == 'pos' else 0)
data3 = data3[~data3.index.isin(['TTN'])]

cpg_list = data3[data3['Is_positive'] == 1].index.values
other_list = data3[data3['Is_positive'] != 1].index.values


# XGboost hyperparameter
xgb_model = XGBClassifier(
    n_estimators=50,         
    max_depth=3,             
    learning_rate=0.01,      
    subsample=0.8,           
    colsample_bytree=0.8,    
    random_state=35       
)

def process_iteration(i):
    DF_gene = []
    DF_score = []
    M2 = []
    KKKK = []
    
    selected_gene_other = np.random.choice(other_list, len(cpg_list), replace=False)
    fin_gene = np.concatenate((cpg_list, selected_gene_other))
    data4 = data3[data3.index.isin(fin_gene)].reset_index()

    X = data4[np.concatenate((['Gene'], x_col))]
    y = data4[['Is_positive']]

    sss = StratifiedShuffleSplit(n_splits=second_repeat, test_size=0.2)
    j = 0

    for train, test in sss.split(X[X.columns[1:]], y['Is_positive']):
        # Dividir datos en entrenamiento y prueba
        X_train_SS = X.iloc[train]
        y_train_SS = y.iloc[train]
        X_test_SS = X.iloc[test]
        y_test_SS = y.iloc[test]
        
        # Entrenar el modelo con los parámetros más usados
        xgb_model.fit(X_train_SS[X_train_SS.columns[1:]], y_train_SS['Is_positive'].values.ravel())
    
        # AUC calculation
        ML_proba = xgb_model.predict_proba(X_test_SS[X_test_SS.columns[1:]])
        malignant_probs = ML_proba[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test_SS, malignant_probs, pos_label=1)
        roc_auc = auc(fpr, tpr)
        auc_df = pd.DataFrame(data={'FPR': fpr, 'TPR': tpr, 'iters': f'{i}_{j}'})
        M2.append(auc_df)      
    
        # PREDICTION and SCORE using best_model
        ML_pred = xgb_model.predict(X_test_SS[X_test_SS.columns[1:]])
        precision = precision_score(y_test_SS, ML_pred)
        recall = recall_score(y_test_SS, ML_pred)
        f1score = f1_score(y_test_SS, ML_pred)
        accuracy = accuracy_score(y_test_SS, ML_pred)
    
        # Tabla de predicciones
        daf = pd.DataFrame(data={
            'Gene': data4[['Gene', 'Is_positive']].iloc[test].values[:, 0],
            'iters': np.repeat(f'{i}_{j}', len(test)),
            'True_label': data4[['Gene', 'Is_positive']].iloc[test].values[:, 1].astype(int),
            'ML_label': ML_pred.astype(int)
        })
    
        # Tabla de scores
        DDD = pd.DataFrame(data={
            'AUC': roc_auc, 'Precision': precision, "Recall": recall, 
            'F1_score': f1score, 'ACCURACY': accuracy
        }, index=[f'{i}_{j}'])
        
        DF_gene.append(daf)
        DF_score.append(DDD)
        j += 1 
        
        KKKK.append([DF_gene, DF_score, M2])
    
    return KKKK

result = Parallel(n_jobs=-1)(delayed(process_iteration)(i) for i in range(first_permu))

re_result_shape = np.array(result, dtype='object').reshape(-1, np.array(result, dtype='object').shape[-1])

GENE = pd.concat([pd.concat(d) for d in re_result_shape[0::4]])
SCORE = pd.concat([pd.concat(d) for d in re_result_shape[1::4]])
AUC_val = pd.concat([pd.concat(d) for d in re_result_shape[2::4]])
GRID_SEARCH_RESULTS = pd.concat([pd.concat(d) for d in re_result_shape[3::4]])

SCORE2 = SCORE #SCORE[(SCORE['AUC']>=0.8)&(SCORE['Precision']>=0.8)]
GENE2 = GENE #GENE[GENE.iters.isin(SCORE2.index)]
AUC_val2 = AUC_val

FIN_GENE = GENE2.groupby('Gene').sum(numeric_only=True)
FIN_GENE['EVENT'] = GENE2.groupby('Gene').count().values[:,1]
FIN_GENE['Likelyhood_CPG'] = FIN_GENE['ML_label'] / GENE2.groupby('Gene').count().values[:,1]


FIN_GENE.to_csv('Result_machine_learning/%s_.ALL.permutation.%s.tsv' %(name, (first_permu*second_repeat)), sep='\t')
SCORE.to_csv('Result_machine_learning/%s.ALL.score.%s.tsv'%(name, (first_permu*second_repeat)) , sep='\t')
AUC_val2.to_csv('Result_machine_learning/%s.ALL.AUC.score.%s.tsv'%(name, (first_permu*second_repeat)) , sep='\t')
GRID_SEARCH_RESULTS.to_csv('Result_machine_learning_gridCV/%s.GRID_SEARCH_RESULTS1.%s.tsv' % (name, (first_permu * second_repeat)), sep='\t')
