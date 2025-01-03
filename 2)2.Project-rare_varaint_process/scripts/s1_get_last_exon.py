import numpy as np
import pandas as pd

vep = pd.read_csv(snakemake.input['TEMPO'], # use file from 008.get_exonic_splicing/out3/s203_VEP_common_TEMPO/VEP.TCGA.{chr}.exonic.filter.cd.seqfiltered.AF.ALL.tempo.KIT.tsv.gz
                    sep='\t', compression='gzip', low_memory=False, 
                     )

clin = 'TCGA-29-2429-10A-01D-1526-09'
vep = vep[vep['SAMPLE']!=clin]

vep['REVEL'] =vep['REVEL'].replace('-',np.nan).apply(pd.to_numeric,errors='coerce')

column_to_AF = ['gnomADe_AF', 'gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_EAS_AF','gnomADe_NFE_AF','gnomADe_SAS_AF','gnomADe_OTH_AF',
           'TCGA_AFR_AF','TCGA_AMR_AF','TCGA_EAS_AF','TCGA_EUR_AF','TCGA_NAN_AF','TCGA_PAS_AF','TCGA_AF']
vep[column_to_AF] = vep[column_to_AF].replace('-',np.nan).apply(pd.to_numeric,errors='coerce')
vep[column_to_AF] = vep[column_to_AF].astype(float)

rare_frequency = float(snakemake.params['maf'])

# 1. get rare variant from all annotation
rare_gnomad_no = vep [    (vep['gnomADe_AF'].isna()) | (vep['gnomADe_AFR_AF'].isna()) | (vep['gnomADe_AMR_AF'].isna()) |
    (vep['gnomADe_EAS_AF'].isna() ) | (vep['gnomADe_NFE_AF'].isna() ) | (vep['gnomADe_OTH_AF'].isna())   ]

rare_gnomad_yes = vep [~vep.ID.isin(rare_gnomad_no.ID)]

rare_gnomad_no = rare_gnomad_no [
    rare_gnomad_no['TCGA_AF'] <0.01
]

rare_gnomad_yes = rare_gnomad_yes[
                (rare_gnomad_yes['gnomADe_AF'] <rare_frequency) &(rare_gnomad_yes['gnomADe_AFR_AF'] <rare_frequency) &
                (rare_gnomad_yes['gnomADe_AMR_AF'] <rare_frequency) &
                (rare_gnomad_yes['gnomADe_EAS_AF'] <rare_frequency) &(rare_gnomad_yes['gnomADe_NFE_AF'] <rare_frequency)
                &(rare_gnomad_yes['gnomADe_OTH_AF'] <rare_frequency)   

]

rare_variant = pd.concat([rare_gnomad_yes, rare_gnomad_no ])




# rare_variant =   vep [
    
#     (vep['gnomADe_AF'] <rare_frequency) &(vep['gnomADe_AFR_AF'] <rare_frequency) &(vep['gnomADe_AMR_AF'] <rare_frequency) &
#     (vep['gnomADe_EAS_AF'] <rare_frequency) &(vep['gnomADe_NFE_AF'] <rare_frequency) &(vep['gnomADe_OTH_AF'] <rare_frequency) &  
    
#     (vep['TCGA_AFR_AF'] <rare_frequency)&(vep['TCGA_AMR_AF'] <rare_frequency)&(vep['TCGA_EAS_AF'] <rare_frequency)&(vep['TCGA_EUR_AF'] <rare_frequency)&
#     (vep['TCGA_NAN_AF'] <rare_frequency)&(vep['TCGA_PAS_AF'] <rare_frequency)&(vep['TCGA_AF'] <rare_frequency)
#               ]




# which included all combination of rare variants -> cutoff of MAF could be alternative
# 

# 2. get clinvar, PTV, Delmis 

clinvar_review = ["practice_guideline",'criteria_provided,_multiple_submitters,_no_conflicts' , 'reviewed_by_expert_panel']


clinvar = rare_variant[(rare_variant['ClinVar_CLNREVSTAT'].isin(clinvar_review)) &
                        (rare_variant['ClinVar_CLNSIG'].isin(['Pathogenic', 'Pathogenic/Likely_pathogenic', 'Likely_pathogenic',
       'Pathogenic/Likely_pathogenic|_risk_factor','Pathogenic|_risk_factor','Pathogenic/Likely_pathogenic|_other',
       'Pathogenic|_drug_response', 'Likely_pathogenic|_drug_response',
       'Pathogenic|_other','Likely_pathogenic|_other','Pathogenic|_drug_response|_other',
        'Likely_pathogenic|_risk_factor','drug_response', 'risk_factor'] )) 

]

ptv = rare_variant[(rare_variant['Consequence'].isin(['frameshift_variant', 'stop_gained', 'start_lost', 'stop_lost',
                                                      'stop_gained,frameshift_variant','stop_gained,inframe_deletion',
                                                      'frameshift_variant,splice_region_variant', 'frameshift_variant,stop_lost,splice_region_variant',
                                                         'start_lost,splice_region_variant', 'frameshift_variant,start_lost,start_retained_variant','inframe_insertion,stop_retained_variant', 'protein_altering_variant'])) |
                                                         (rare_variant['Consequence'].apply(lambda x: 'splice' in x))
]


revel = rare_variant[(rare_variant['REVEL'] > 0.5) & (rare_variant['Consequence'].apply(lambda x: 'missense_variant' in x))]

benign = rare_variant[(rare_variant['ClinVar_CLNREVSTAT'].isin(clinvar_review )) &
    (rare_variant['ClinVar_CLNSIG'].isin(['Likely_benign','Benign/Likely_benign', 'Benign']))]


rare_variant['is_delmis'] = rare_variant['ID'].isin(revel['ID'].tolist())
rare_variant['is_clin'] = rare_variant['ID'].isin(clinvar['ID'].tolist())
rare_variant['is_ptv'] = rare_variant['ID'].isin(ptv['ID'].tolist())
rare_variant['is_benign'] = rare_variant['ID'].isin(benign['ID'].tolist())

rdgv =  rare_variant[(rare_variant['is_clin'] == True) | ((rare_variant['is_benign'] != True) & (rare_variant['is_ptv'] == True)) | ((rare_variant['is_benign'] != True) & (rare_variant['is_delmis'] == True)) ]

rare_variant['is_rdgv'] = rare_variant['ID'].isin(rdgv['ID'].tolist())


#########################################
# check SpliceAI

PTV_other = rare_variant[
    (rare_variant['is_rdgv']==True) & (rare_variant['is_ptv']==True) & (rare_variant['Consequence'].isin(['frameshift_variant', 'stop_gained', 'start_lost', 'stop_lost',
                                                                                                          'stop_gained,frameshift_variant','stop_gained,inframe_deletion','frameshift_variant,splice_region_variant', 'frameshift_variant,stop_lost,splice_region_variant',
                                                         'start_lost,splice_region_variant', 'frameshift_variant,start_lost,start_retained_variant','inframe_insertion,stop_retained_variant', 'protein_altering_variant']))
                         ]
# 'stop_gained,frameshift_variant','stop_gained,inframe_deletion', -> need to include!


PTV_splicing = rare_variant[
    (rare_variant['is_rdgv']==True) & (rare_variant['is_ptv']==True) & 
    (rare_variant['Consequence'].apply(lambda x: 'splice' in x)) 
                            ]

PTV_splicing_acc_don_true = PTV_splicing[PTV_splicing['Consequence'].apply(lambda x: 'splice_donor_variant' in x) |PTV_splicing['Consequence'].apply(lambda x: 'splice_acceptor_variant' in x) ]

PTV_spliceAI = PTV_splicing[~PTV_splicing.ID.isin(PTV_splicing_acc_don_true.ID)]

Splice_above_8 = PTV_spliceAI[PTV_spliceAI['SpliceAI_pred'].apply(lambda x: (x != '-' and len(x.split('|')) >= 5 and  any(float(x.split('|')[i]) >= 0.8 for i in range(1, 5))
                                                                             ))]

# PTV merge after measuring SPliceAI
PTV_assemble = pd.concat([PTV_other, PTV_splicing_acc_don_true, Splice_above_8])

PTV_assemble_exonic  = PTV_assemble[PTV_assemble['EXON'].apply(lambda x: "-" not in x)]
PTV_assemble_intronic = PTV_assemble[PTV_assemble['INTRON'].apply(lambda x: "-" not in x)]

# if last exon does not exist?
# 

# check last exon in exon
PTV_assemble_exonic_yes = PTV_assemble_exonic[PTV_assemble_exonic['EXON'].apply(lambda x: ( x.split('/')[0] == x.split('/')[1] )) ]
PTV_assemble_exonic_no = PTV_assemble_exonic[PTV_assemble_exonic['EXON'].apply(lambda x: (x.split('/')[0] != x.split('/')[1] ) )]

PTV_assemble_exonic_yes['is_last_exon_intron'] = True
PTV_assemble_exonic_no['is_last_exon_intron'] = False

# check last exon in intronic
PTV_assemble_intronic_yes = PTV_assemble_intronic[PTV_assemble_intronic['INTRON'].apply(lambda x: ( x.split('/')[0] == x.split('/')[1] ))  ]
PTV_assemble_intronic_no = PTV_assemble_intronic[PTV_assemble_intronic['INTRON'].apply(lambda x: (x.split('/')[0] != x.split('/')[1] ) )]

PTV_assemble_intronic_yes['is_last_exon_intron'] = True
PTV_assemble_intronic_no['is_last_exon_intron'] = False




PTV_last_fin = pd.concat([PTV_assemble_exonic_yes,PTV_assemble_exonic_no,PTV_assemble_intronic_yes,PTV_assemble_intronic_no])

PTV_last_fin.to_csv(snakemake.output['OUT'], sep='\t',compression='gzip', index=False)

rare_variant.to_csv(snakemake.output['ALL'],sep='\t', compression='gzip',index=False)