import numpy as np
import pandas as pd

ptex_combined = pd.read_csv(snakemake.input['ptex_combined'],
                            compression='gzip', sep='\t', low_memory=False)

tcga = ['ACC','BLCA','BRCA','CESC','CHOL','COADREAD','DLBC','ESCA','GBM','HNSC','KICH','KIRP','KIRC',
'LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','SARC','SKCM','STAD','THCA','THYM','TGCT','UCEC','UCS','UVM']

gtex = ['AdrenalGland','Bladder','Breast_MammaryTissue','Vagina','mean_proportion','Colon_Transverse','WholeBlood','Esophagus_Mucosa',
'Brain_Substantianigra','mean_proportion','Kidney_Cortex','Kidney_Cortex','Kidney_Cortex','WholeBlood','Brain_Amygdala',
'Liver','Lung','Lung','mean_proportion','Ovary','Pancreas','Pituitary','Prostate','mean_proportion','Skin_SunExposed_Lowerleg_',
'Stomach','Thyroid','Testis','mean_proportion','FallopianTube','Cervix_Endocervix','mean_proportion']

tcga_gtex = dict(map(lambda i,j : (i,j) , tcga,gtex))

mm = ptex_combined

df_s_column = ['ID', 'POS', 'FILTER','AC','gnomAD_FILTER','SAMPLE','GT','AD','RD','DP','GQ','VAF','CANCER','ANCESTRY',
               'CENTER', 'GENDER', 'SEQ_KIT', 'SAM_LOWQUAL', 'BI_AL','WUSM_AL','BCM_AL','WTSI_AL','CENTER_AL','BCM_AF','BI_AF','WTSI_AF','WUSM_AF','CENTER_AF',       
               'AFR_AL', 'AMR_AL', 'EAS_AL', 'EUR_AL', 'NAN_AL', 'PAS_AL', 'TCGA_AL','TCGA_AFR_AF', 'TCGA_AMR_AF','TCGA_EAS_AF', 'TCGA_EUR_AF', 'TCGA_NAN_AF', 
               'TCGA_PAS_AF', 'TCGA_AF', 'BCM_KIT_AL','BI_KIT_AL', 'SeqCap2_AL', 'SeqCap3_AL','WTSI_KIT_AL','WUSM_KIT38_AL', 'WUSM_KIT18_AL','KIT_AL',
                'BCM_KIT_AF','BI_KIT_AF','SeqCap2_AF' ,'SeqCap3_AF','WTSI_KIT_AF', 'WUSM_KIT38_AF','WUSM_KIT18_AF', 'KIT_AF',  
                'Gene','Feature','Feature_type','Consequence','SYMBOL','SYMBOL_SOURCE','BIOTYPE','CANONICAL','SIFT','PolyPhen','EXON','INTRON',  
                'gnomADe_AF','gnomADe_AFR_AF', 'gnomADe_AMR_AF',  'gnomADe_ASJ_AF', 'gnomADe_EAS_AF',  'gnomADe_FIN_AF', 
                'gnomADe_NFE_AF','gnomADe_OTH_AF','gnomADe_SAS_AF','SpliceAI_pred','REVEL','ClinVar', 'ClinVar_CLNSIG','ClinVar_CLNREVSTAT', 
                'ClinVar_CLNDN','is_delmis','is_clin', 'is_ptv',  'is_benign','is_rdgv','ensg', 'symbol',]

new_column  = ['ID', 'POS', 'FILTER','AC','gnomAD_FILTER','SAMPLE','GT','AD','RD','DP','GQ','VAF','CANCER','ANCESTRY',
               'CENTER', 'GENDER', 'SEQ_KIT', 'SAM_LOWQUAL', 'BI_AL','WUSM_AL','BCM_AL','WTSI_AL','CENTER_AL','BCM_AF','BI_AF','WTSI_AF','WUSM_AF','CENTER_AF',       
               'AFR_AL', 'AMR_AL', 'EAS_AL', 'EUR_AL', 'NAN_AL', 'PAS_AL', 'TCGA_AL','TCGA_AFR_AF', 'TCGA_AMR_AF','TCGA_EAS_AF', 'TCGA_EUR_AF', 'TCGA_NAN_AF', 
               'TCGA_PAS_AF', 'TCGA_AF', 'BCM_KIT_AL','BI_KIT_AL', 'SeqCap2_AL', 'SeqCap3_AL','WTSI_KIT_AL','WUSM_KIT38_AL', 'WUSM_KIT18_AL','KIT_AL',
                'BCM_KIT_AF','BI_KIT_AF','SeqCap2_AF' ,'SeqCap3_AF','WTSI_KIT_AF', 'WUSM_KIT38_AF','WUSM_KIT18_AF', 'KIT_AF',  
                'Gene','Feature','Feature_type','Consequence','SYMBOL','SYMBOL_SOURCE','BIOTYPE','CANONICAL','SIFT','PolyPhen','EXON','INTRON',  
                'gnomADe_AF','gnomADe_AFR_AF', 'gnomADe_AMR_AF',  'gnomADe_ASJ_AF', 'gnomADe_EAS_AF',  'gnomADe_FIN_AF', 
                'gnomADe_NFE_AF','gnomADe_OTH_AF','gnomADe_SAS_AF','SpliceAI_pred','REVEL','ClinVar', 'ClinVar_CLNSIG','ClinVar_CLNREVSTAT', 
                'ClinVar_CLNDN','is_delmis','is_clin', 'is_ptv',  'is_benign','is_rdgv','ensg','symbol', "PTEX"]



big=[]
for index, row in mm.iterrows():
    small=[]
    #print(index)
    keyss = mm.loc[index,'CANCER']
    valuess = mm.loc[index,tcga_gtex.get(keyss)]
    for k in df_s_column:
        small.append(mm.loc[index,k])
    #print(valuess)
    #print(small)
    small.append(valuess)
    #small.append(mm.loc[index,['ID','SAMPLE','CANCER','ANCESTRY','GENDER','CENTER','Gene.ensGene']].values)
    #small.append(valuess)
    #break
    big.append(small)

FINAL = pd.DataFrame(data=big, columns=new_column)


FINAL2 = FINAL
FINAL2  = FINAL2.astype({'PTEX':'float'})

FINAL2['is_ptex'] = FINAL2['PTEX'].apply(lambda x: x > 0.1 )

FINAL2.to_csv(snakemake.output['merged'],
              sep='\t', compression='gzip', index=False)