import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, precision_score, recall_score, f1_score, accuracy_score
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

# Cargar los datos
data = pd.read_csv('./TCGA.integrative.model.input.by.rdgv.noSOD12.with.new_PPI', sep='\t', index_col=0)

# Rellenar valores faltantes
data['avg_CPG'] = data['avg_CPG'].fillna(np.max(data['avg_CPG']))
data['avg_DNA_repair'] = data['avg_DNA_repair'].fillna(np.max(data['avg_DNA_repair']))
data['avg_PIK_mTOR'] = data['avg_PIK_mTOR'].fillna(np.max(data['avg_PIK_mTOR']))
data['avg_Cell_Cycle'] = data['avg_Cell_Cycle'].fillna(np.max(data['avg_Cell_Cycle']))
data['avg_ALL'] = data['avg_ALL'].fillna(np.max(data['avg_ALL']))
data['distanceAll_path'] = data['distanceAll_path'].fillna(np.max(data['distanceAll_path']))
data = data.fillna(0)



# Normalizar solo las columnas de características
features = ['case_con_zvalues', 'two_hit_zvalues', 'RDGV_by_nonRDGV', 'Tau_score',
       'oe_lof_upper', 'PTM_pvalue', 'PTV_frequency', 'amplification.freq',
       'deletion.freq', 'clin_zvalue', 'distanceAll_path',
'Bidirectional_Deletion', 'Edition_efficiency', 'Guide_Abundance',
       'Insertion', 'Insertion_Deletion', 'PAM_Proximal_Deletion', 'Is_positive']


scaler = StandardScaler()
data[features] = scaler.fit_transform(data[features])

# Convertir la columna de etiquetas a valores binarios
data['Is_positive'] = data['Is_positive'].apply(lambda x: 1 if x == 'pos' else 0)

# Remover el índice 'TTN' si está presente
data = data[~data.index.isin(['TTN'])]

# Definir listas de índices de cada clase
cpg_list = data[data['Is_positive'] == 1].index.values
other_list = data[data['Is_positive'] != 1].index.values

# Configurar regresión logística y GridSearchCV para encontrar los mejores parámetros
logreg = LogisticRegression(max_iter=1000)


# logistic regression hyperparameter
param_grid = {
    'C': np.logspace(-3, 3, 7),
    'solver': ['liblinear', 'saga'],
    'penalty': ['l1', 'l2', 'elasticnet'],
    'l1_ratio': [0.1, 0.5, 0.7, 0.9]
}


first_permu = 2000 
second_repeat = 5
x_col = features

def process_iteration(i):
    DF_gene = []
    DF_score = []
    M2 = []
    GS_results = []

    selected_gene_other = np.random.choice(other_list, len(cpg_list), replace=False)
    fin_gene = np.concatenate((cpg_list, selected_gene_other))
    data4 = data.loc[fin_gene].reset_index()
    
    X = data4[x_col]
    y = data4['Is_positive']
    
    sss = StratifiedShuffleSplit(n_splits=second_repeat, test_size=0.2)
    j=0

    for train, test in sss.split(X, y):
        X_train, X_test = X.iloc[train], X.iloc[test]
        y_train, y_test = y.iloc[train], y.iloc[test]
        
        # Aplicar GridSearchCV para optimizar la regresión logística
        grid_search = GridSearchCV(logreg, param_grid, cv=5, scoring='accuracy', n_jobs=-1, verbose=0)
        grid_search.fit(X_train, y_train)
        best_model = grid_search.best_estimator_

        GS_results.append(pd.DataFrame([grid_search.best_params_], index=[f'{i}_{j}']))

        # AUC
        y_pred_proba = best_model.predict_proba(X_test)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
        roc_auc = auc(fpr, tpr)
        auc_df = pd.DataFrame({'FPR': fpr, 'TPR': tpr, 'iters': f'{i}_{j}'})
        M2.append(auc_df)
        
        # Predicción y métricas de desempeño
        y_pred = best_model.predict(X_test)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        f1score = f1_score(y_test, y_pred)
        accuracy = accuracy_score(y_test, y_pred)

        # Tabla de probabilidades
        daf = pd.DataFrame({
            'Gene': data4['Gene'].iloc[test],
            'iters': f'{i}_{j}',
            'True_label': y_test.values,
            'ML_label': y_pred
        })
        
        # Tabla de puntuaciones
        DDD = pd.DataFrame({
            'AUC': roc_auc, 'Precision': precision, "Recall": recall, 
            'F1_score': f1score, 'ACCURACY': accuracy
        }, index=[f'{i}_{j}'])
        
        DF_gene.append(daf)
        DF_score.append(DDD)
        j+=1
        
    return DF_gene, DF_score, M2, GS_results

# Ejecutar las permutaciones en paralelo
result = Parallel(n_jobs=-1)(delayed(process_iteration)(i) for i in range(first_permu))

re_result_shape =np.array(result,dtype='object').reshape(-1, np.array(result,dtype='object').shape[-1])
GENE = pd.concat([pd.concat(d) for d in re_result_shape[0::4]])
SCORE = pd.concat([pd.concat(d) for d in re_result_shape[1::4]])
AUC_val = pd.concat([pd.concat(d) for d in re_result_shape[2::4]])
GRID_SEARCH_RESULTS = pd.concat([pd.concat(d) for d in re_result_shape[3::4]])



# Guardar los resultados
name = 'LogReg_impu'
FIN_GENE = GENE.groupby('Gene').sum(numeric_only=True)
FIN_GENE['EVENT'] = GENE.groupby('Gene').count().values[:, 1]
FIN_GENE['Likelyhood_CPG'] = FIN_GENE['ML_label'] / GENE.groupby('Gene').count().values[:, 1]

FIN_GENE.to_csv('Result_machine_learning/%s_.permutation1.%s.tsv' %(name, (first_permu*second_repeat)), sep='\t')
SCORE.to_csv('Result_machine_learning/%s.score1.%s.tsv'%(name, (first_permu*second_repeat)) , sep='\t')
AUC_val.to_csv('Result_machine_learning/%s.AUC.score1.%s.tsv'%(name, (first_permu*second_repeat)) , sep='\t')
GRID_SEARCH_RESULTS.to_csv('Result_machine_learning_gridCV/%s.GRID_SEARCH_RESULTS1.%s.tsv' % (name, (first_permu * second_repeat)), sep='\t', index=False)



# annex to draw AUC figure, but not usable for this step
# max_auc_iter = AUC_val.groupby('iters').apply(lambda x: auc(x['FPR'], x['TPR'])).idxmax()
# best_auc_df = AUC_val[AUC_val['iters'] == max_auc_iter]

# # Graficar la curva ROC para el AUC más alto
# plt.figure(figsize=(10, 6))
# plt.plot(best_auc_df['FPR'], best_auc_df['TPR'], color='blue', label=f'Best AUC (iteration {max_auc_iter})')
# plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Best ROC Curve for Logistic Regression')
# plt.legend(loc="lower right")
# plt.savefig('/storage/scratch01/groups/gcc/alex_pm_project_24/005.MACHINE_LEARNING/003.ML_tries/best_roc_curve.png')

