import numpy as np

import pandas as pd

# ptex problem solution

ptex = pd.read_csv(snakemake.input['ptex'],
                   compression='gzip', sep='\t', low_memory=False)

delmis = pd.read_csv(snakemake.input['delmis'],
                    compression='gzip', sep='\t', low_memory=False)


delmis['CHR'] = delmis['ID'].apply(lambda x: x.split('_')[0]) 

#ptex

ptex['chr'], ptex['pos'] = ptex['locus'].apply(lambda x: x.split(':')[0]), ptex['locus'].apply(lambda x: x.split(':')[1])
ptex['pos'] = ptex['pos'].astype(int)


#ptex = ptex.drop(columns=['locus'])

ptex2 = ptex[['chr', 'pos', 'ensg', 'symbol', 'Cells_Transformedfibroblasts', 'Prostate', 'Spleen',
                'Brain_FrontalCortex_BA9_', 'SmallIntestine_TerminalIleum',
                'MinorSalivaryGland', 'Artery_Coronary', 'Skin_SunExposed_Lowerleg_',
                'Cells_EBV_transformedlymphocytes', 'Brain_Hippocampus',
                'Esophagus_Muscularis', 'Brain_Nucleusaccumbens_basalganglia_',
                'Artery_Tibial', 'Brain_Hypothalamus', 'Adipose_Visceral_Omentum_',
                'Cervix_Ectocervix', 'Brain_Spinalcord_cervicalc_1_',
                'Brain_CerebellarHemisphere', 'Nerve_Tibial', 'Breast_MammaryTissue',
                'Liver', 'Skin_NotSunExposed_Suprapubic_', 'AdrenalGland', 'Vagina',
                'Pancreas', 'Lung', 'FallopianTube', 'Pituitary', 'Muscle_Skeletal',
                'Colon_Transverse', 'Artery_Aorta', 'Heart_AtrialAppendage',
                'Adipose_Subcutaneous', 'Esophagus_Mucosa', 'Heart_LeftVentricle',
                'Brain_Cerebellum', 'Brain_Cortex', 'Thyroid', 'Brain_Substantianigra',
                'Kidney_Cortex', 'Uterus', 'Stomach', 'WholeBlood', 'Bladder',
                'Brain_Anteriorcingulatecortex_BA24_', 'Brain_Putamen_basalganglia_',
                'Brain_Caudate_basalganglia_', 'Colon_Sigmoid', 'Cervix_Endocervix',
                'Ovary', 'Esophagus_GastroesophagealJunction', 'Testis',
                'Brain_Amygdala', 'mean_proportion']]

ptex_delmis = delmis.merge(ptex2, left_on=['CHR','POS','Gene','SYMBOL'], right_on =['chr','pos','ensg','symbol'], how='left')

#del(ptex)
#del(delmis)

ptex_delmis.to_csv(snakemake.output['ptex_delmis'],
                   sep='\t', compression='gzip', index=False)

