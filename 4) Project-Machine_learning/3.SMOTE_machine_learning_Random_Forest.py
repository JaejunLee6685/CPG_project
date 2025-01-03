import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.ensemble import RandomForestClassifier  # Cambiado a RandomForest
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from joblib import Parallel, delayed

# too old to update, if SMOTE make sense, it is better to update with new version

# Cargar los datos
data = pd.read_csv('./TCGA.integrative.model.input.by.rdgv.noSOD12.with.new_PPI.tsv', sep='\t')

# Reemplazar valores nulos por el máximo de cada columna
data['avg_CPG'] = data['avg_CPG'].fillna(np.max(data['avg_CPG']))
data['avg_DNA_repair'] = data['avg_DNA_repair'].fillna(np.max(data['avg_DNA_repair']))
data['avg_PIK_mTOR'] = data['avg_PIK_mTOR'].fillna(np.max(data['avg_PIK_mTOR']))
data['avg_Cell_Cycle'] = data['avg_Cell_Cycle'].fillna(np.max(data['avg_Cell_Cycle']))
data['avg_ALL'] = data['avg_ALL'].fillna(np.max(data['avg_ALL']))
data = data.fillna(0)

# Listas de genes en clases positivas y negativas
cpg_list = data[data['Is_positive'] == 'pos']['Gene'].values  # Clase minoritaria
other_list = data[data['Is_positive'] == 'neg']['Gene'].values  # Clase mayoritaria

repeat = 2000
splits = 5
INDEXXX = data.columns[1:-2]  # Características

np.random.seed(35)

# Inicialización del clasificador con RandomForest
classifier = RandomForestClassifier(n_estimators=1000, random_state=35) # default random forest

def process_iteration(i):
    # No submuestreamos la clase mayoritaria. Usamos todos los ejemplos de other_list (15,000)
    fin_gene = np.concatenate((cpg_list, other_list))  # Combinamos todos los genes
    
    data3 = data[data.Gene.isin(fin_gene)]  # Conjunto de datos combinado

    # Separar características (X) y variable objetivo (y)
    X = data3[data.columns[:-2]]
    y = data3[['Is_positive']]

    # Normalizar características numéricas
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X[X.columns[1:]])  # Normalizar sin incluir la columna 'Gene'

    # Stratified Shuffle Split para mantener las proporciones de clases
    sss = StratifiedShuffleSplit(n_splits=splits, test_size=0.2, random_state=35)
    KKKK = []
    
    for train_idx, test_idx in sss.split(X_scaled, y['Is_positive']):
        X_train, y_train = X_scaled[train_idx], y.iloc[train_idx]

        # Aplicar SMOTE solo en los datos de entrenamiento (sobre la clase minoritaria cpg_list)
        smote = SMOTE(random_state=35) # seed in SMOTE function
        X_train_SMOTE, y_train_SMOTE = smote.fit_resample(X_train, y_train['Is_positive'])

        # Entrenar el clasificador con los datos balanceados
        sel = SelectFromModel(classifier)
        sel.fit(X_train_SMOTE, y_train_SMOTE.values.ravel())

        # Extraer la importancia de las características desde el modelo entrenado en SelectFromModel
        importances = sel.estimator_.feature_importances_

        # Verificar qué características fueron seleccionadas (True/False)
        selected_features_mask = sel.get_support()
        t_or_f = selected_features_mask.astype(int)
        
        KKKK.append(importances)
        KKKK.append(t_or_f)
    
    return KKKK

# Ejecutar el bucle en paralelo
results = Parallel(n_jobs=-1)(delayed(process_iteration)(i) for i in range(repeat)) # n=-1 is not recommended

# Reshape de los resultados
reshaped_array = np.array(results).reshape(-1, np.array(results).shape[-1])
even_indices = reshaped_array[0::2]  # Importancias
odd_indices = reshaped_array[1::2]   # Selección de características

# Guardar los resultados en archivos
even_df = pd.DataFrame(data=even_indices, columns=INDEXXX)
odd_df = pd.DataFrame(data=odd_indices, columns=INDEXXX)

even_df.to_csv('./Result_feature_importance/RF_impu.importance.SMOTE.tsv.gz', sep='\t', compression='gzip')
odd_df.to_csv('./Result_feature_importanceRF_impu.small_event.SMOTE.tsv.gz', sep='\t', compression='gzip')
