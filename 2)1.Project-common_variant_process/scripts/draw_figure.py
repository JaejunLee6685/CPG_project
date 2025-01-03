import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

# import basic file and preprocessing step
clin = pd.read_csv(snakemake.input['clin'],sep='\t')

pca = pd.read_csv(snakemake.input['eigenvec'], low_memory=False, sep=' ', header=None)
pca.columns = ['PC'+str(x-1) for x in pca.columns]
pca = pca.drop(columns='PC-1')
pca = pca.rename(columns={'PC0':'SAMPLE'})

# use as the sample low quality stat
merged_pca = pca.merge(clin, left_on=['SAMPLE'],right_on=['SAMPLE'],).sort_values('SAM_LOWQUAL', ascending=True)
#merged_pca2 = merged_pca.sort_values('SAM_LOWQUAL', ascending=True)

eigen_val = pd.read_csv(snakemake.input['eigenval'],sep=' ', header=None)
cova = pd.read_csv(snakemake.input['eigencov'],sep=' ', header=None)
pve = pd.DataFrame(data=(eigen_val/np.sum(cova)*100).values, index=merged_pca.columns[1:len(merged_pca.columns)-9], 
                   columns=['PC_percent'])

# shape of graphs
def make_alpha(x):
    return to_rgba(x,1)
def make_shape(x):
    return 'o'
def take_alpha(x):
    return x[3]

# color by manual
pops = dict(zip(['EUR', 'EAS', 'AFR' ,'NAN', 'AMR', 'PAS'], ["#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A","#666666", "#E6AB02", "#A6761D"]) )
centre = dict(zip(['BI', 'WUSM' ,'BCM' ,'WTSI'],("#66C2A5", "#FC8D62" ,"#E78AC3","#8DA0CB") ))
origin = dict(zip( ['BLOOD' ,'NORMAL' ,'TUMOR' ,'BUCCAL', 'BONE'] ,("#E41A1C", "#A65628", "#984EA3", '#377EB8', '#000000')))
kits = dict(zip(['BI_refseq_plus_3', 'WUSM_38', 'WUSM_2_1', 'WUSM_50', 'BI_other' ,'BI_tcga_6k', 'SeqCap_V3', 'WUSM_hg18' ,'BCM_Nimble_VCRome', 'Unknown' ,'WTSI_exon_V5', 'SeqCap_V2'] ,  
                ("#F07E4D","#FD6F88","#E6AB02","#82AD00","#00C1A7","#619CFF","#FF61CA","#00B5ED","#EF67EB","#8DA0CB","#CE9500","#00AEF9" )))
cancers = dict(zip( ['GBM' ,'OV', 'LUAD' ,'LUSC' ,'PRAD' ,'UCEC' ,'BLCA' ,'TGCT' ,'ESCA', 'PAAD', 'KIRP',
 'LIHC' ,'CESC','SARC' ,'BRCA', 'THYM' ,'MESO', 'COADREAD', 'STAD', 'SKCM' ,'CHOL',
 'KIRC', 'THCA', 'HNSC' ,'LAML', 'LGG', 'DLBC', 'KICH', 'UCS', 'ACC' ,'PCPG' ,'UVM'],
        ("#F8766D" ,"#F07E4D", "#E78619" ,"#DB8E00" ,"#CE9500" ,"#BF9C00" ,"#AEA200",
                "#9AA800" ,"#82AD00" ,"#64B200" ,"#32B600" ,"#00BA38", "#00BD5C" ,"#00BF78",
                "#00C091" ,"#00C1A7" ,"#00C0BB" ,"#00BECD", "#00BADE" ,"#00B5ED", "#00AEF9",
                "#00A6FF" ,"#619CFF" ,"#9191FF", "#B385FF", "#CD79FF", "#E16FF8", "#EF67EB",
                "#F962DB" ,"#FF61CA" ,"#FF63B6" ,"#FF68A0" ,"#FD6F88")      ))
genders = dict(zip(['M', 'F'],('#FF007F','#4B0082') ))
seq_low =dict(zip( [False , True] ,('#898989', '#FF0000') ))
sam_low =dict(zip( [False , True] ,('#898989', '#FF0000') ))

legend_name = ['ANCESTRY','CENTER','ORIGIN','SEQ_KIT','CANCER','GENDER','SAM_LOWQUAL','SEQ_LOWQUAL']


col_list_dict =[]
for i in [pops,centre,origin,kits,cancers,genders,sam_low,seq_low]:
    col_list_dict.append(dict(zip(i.keys(), map(make_alpha,i.values()))))

legend_color = []
for i in [pops,centre,origin,kits,cancers,genders,sam_low,seq_low]:
    legend_color.append([k for k in i.values()])

legend_alpha = []
for i in col_list_dict:
    legend_alpha.append([k for k in map(take_alpha,i.values() ) ])

PCs= np.array([[['PC1','PC2'],['PC1','PC3'],['PC1','PC4']],
              [['PC2','PC3'],['PC2','PC4'],['PC3','PC4']]])
HUES = np.array([[['PC1','PC2'],['PC1','PC3'],['PC1','PC4']],
              [['PC2','PC3'],['PC2','PC4'],['PC3','PC4']]])

LABELS = []
for i in col_list_dict:
    LABELS.append([k for k in i.keys()])

MARKERS = []
for i in col_list_dict:
    MARKERS.append(dict(zip(i.keys(), map(make_shape,i.values() ) ) ) ) 



# ORIGINAL plot 
for hues,pals,king,qu,ss,mm,FFF in zip(legend_name,col_list_dict,legend_color,legend_alpha,LABELS,MARKERS,snakemake.output['commons'][0:] ) :
    TOKEN=1
    position_bbox = (1,0.6)
    fontsize='20'
    markersize=20
    title_fontsize='20'
    if hues == "CANCER":
        TOKEN =1
        position_bbox = (1.,1)
    elif hues =='SEQ_KIT':
        TOKEN =1
        position_bbox = (1,0.6)
        fontsize='15'
        markersize=15
        title_fontsize='15'
    else:
        pass
    
    fig, axs = plt.subplots(ncols=3,nrows=2, figsize=(30,15))
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    for k in [0,1]:
        for i in [0,1,2]:
            #print(i)
            sns.scatterplot(data=merged_pca, x=PCs[k][i][0],y=PCs[k][i][1], ax=axs[k,i],
                       hue=hues, palette=pals, legend=True, style=hues, markers=mm)
            axs[k,i].set_xlabel(PCs[k][i][0]+' ('+ str(np.around(pve.loc[PCs[k][i][0]][0],2) )+'%' +')' ,fontsize=20)
            axs[k,i].set_ylabel(PCs[k][i][1]+ ' ('+ str(np.around(pve.loc[PCs[k][i][1]][0],2)) +'%' +')',fontsize=20)
            axs[k,i].tick_params(axis='x',size=5,labelsize=15)
            axs[k,i].tick_params(axis='y',size=5,labelsize=15)
            axs[k,i].get_legend().remove()
    h, l = axs[k,i].get_legend_handles_labels()   
    
    LEG=fig.legend(h,l,title=hues,bbox_to_anchor=position_bbox,fontsize=fontsize,ncol=TOKEN,title_fontsize=title_fontsize,)

    for M in range(len(LEG.legend_handles)):
        #LEG.legendHandles[M].set_color(king[M])
        LEG.legend_handles[M].set(alpha = qu[M], markersize=markersize, )
        #LEG.legend_handles[M].set_sizes([500])
    fig.suptitle('\nTCGA PCA analysis %s' %hues, fontsize=35)
    

    plt.savefig(FFF, )