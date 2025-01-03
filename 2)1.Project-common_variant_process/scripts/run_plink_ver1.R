library('plinkQC')
library('ggplot2')

#setwd('/storage/scratch01/groups/gcc/LJJ_ML_project_2024/006.TCGA_VCF_FILTRATION/009.PCA_common') # necessary?

vcf_dir = snakemake@input['vcf'] # dir of vcf file
#folds = snakemake@params['folds']
dir_folds = snakemake@params['dir_folds']
# make directory by snakemake.parameters 

# make "out" dir
dir.create(file.path(dir_folds), recursive=T, showWarnings = F)

# make dir out/snakemake@params['dir_folds']

##################
plink_dir = 'plink' # find out PLINK directory for sure
##################


plink_prefix = dir_folds

#sprintf('s_fold_file_%s', folds)
 # how to make parameter like fold or directory

ld_dir = '001.data/high-LD-regions-hg19-GRCh37.bed'
clin_data = '001.data/TCGA_SAMPLE_CLIN.txt'
ld_par = '50 5 0.1'
hwe = '1e-6'
file_pre = 'TCGA_fold'
cluster = 'ANCESTRY'
outlier = '/001.data/out.tsv'

#outlier = './003.result/ENSG_all_common/NON_EUR_SAMPLE.tsv'

#################
## 1. make BED ##
# without outlier
#system(command = sprintf('%s --vcf %s --make-bed --out %s/1_1_init',plink_dir, vcf_dir, plink_prefix))

# original code
system(command = sprintf('%s --vcf %s --make-bed --out %s/1_init',plink_dir, vcf_dir, plink_prefix))

## 1-1. remove outlier
system(command = sprintf('%s --bfile %s/1_init --remove %s --make-bed --out %s/1_1_init', plink_dir, plink_prefix, outlier, plink_prefix))

#################

#################
## 2. Heterozygosity Rate Calling ##

fail_het_imiss <- check_het_and_miss(indir = plink_prefix,
                                     name = '1_1_init',
                                    imissTh = 0.01, # 0.03
                                     hetTh = 3,
                                     interactive = F,
                                     label_fail = F, path2plink = plink_dir,
                                     showPlinkOutput = F)

write.table(data.frame(fail_het_imiss$fail_het$FID, fail_het_imiss$fail_het$IID),
            file = paste0(plink_prefix, '/2_heterotest.fail.IDs'), sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

system(command = sprintf('%s --bfile %s/1_1_init --remove %s/2_heterotest.fail.IDs --make-bed --out %s/2_P_hetero',
                         plink_dir, plink_prefix, plink_prefix, plink_prefix))
#################

#################
## 3. LD computation and high LD region ##
system(command=sprintf('%s --bfile %s/2_P_hetero --exclude %s --range --indep-pairwise %s --make-bed --out %s/3_P_ld',
                       plink_dir, plink_prefix, ld_dir, ld_par, plink_prefix))
#################

#################
## 4. Autosome and HWE equilibrium ##
system(command = sprintf('%s --bfile %s/3_P_ld --extract %s/3_P_ld.prune.in --autosome --hwe %s --make-bed --out %s/4_hweq',
                         plink_dir, plink_prefix, plink_prefix, hwe, plink_prefix ))
#################

#################
## 5. PCA  ##
system(command = sprintf('%s --bfile %s/4_hweq --make-rel square0 --pca 8000 --out %s/5_pca',
                         plink_dir, plink_prefix, plink_prefix))

system(command = sprintf("awk '{print $NR}' %s/5_pca.rel > %s/5_pca.diagonal.covariance", plink_prefix, plink_prefix))



# under the code represents the result figure
# but, it is not well constructed code
# we drew the figure with python not R.
# #################################################
# pca <- read.table(paste0(plink_prefix,'/5_pca.eigenvec'), header =FALSE )
# pca <- pca[,-2]
# names(pca)[1] <- "TCGA_PatientID"
# names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# pop_info <- read.table(clin_data, header = 1)
# pop_info <-pop_info[,-2]
# names(pop_info)[1] <- "TCGA_PatientID"

# pca <- merge(pca, pop_info, by = 'TCGA_PatientID') 

# ##############################
# eigenval <- scan(paste0(plink_prefix,'/5_pca.eigenval'))
# covarmatrix_diag <- scan(paste0(plink_prefix,'/5_pca.diagonal.covariance'))

# pve <- data.frame(PC = 1:length(eigenval), p.v.e = eigenval/sum(covarmatrix_diag)*100)

# pve_sub <- pve[1:20, ] 

# pve_plot <- ggplot(pve_sub, aes(PC, p.v.e)) + geom_bar(stat = 'identity', fill = '#7E52A0') + 
#   geom_text(aes(label = round(p.v.e, 2)), position = position_dodge(.9), vjust = -0.5, size = 2.8) +
#   ylab("Percentage variance explained") + ggtitle("Percentage Variance Explained of each Principal Component") + 
#   theme_update(plot.title = element_text(hjust = 0.5)) +
#   theme(axis.title.y =  element_text(size = 14, vjust = 1.5)) +
#   scale_x_continuous(breaks = 1:nrow(pve_sub))

# ggsave(plot = pve_plot, path = paste0(plink_prefix , '/'), height = 5, width = 8,
#        filename = paste0(file_pre, '_population_pve.png'), device = png)


# # how many pc will be drawn
# #components <- combn(seq(1, 4), m = 2) # PC 1, PC 2 

# a <- 1 #components[, comb][1]
# b <- 1 #components[, comb][2]
# j<-0
# for (a in c(1,1,1,2,2,3) ){
#   j =j+1
#     if (a==2){
#     b=2

#   }
#   else if (a==3){
#     b=3
#   }
#   if (j ==5) {
#     b=3
#   }
#   b=b+1

# color <- 'Dark2'
# legend <- 'Population'
# plot_subtitle = ''
# plot_title = 'Common Variants (MAF >= 5%) Across all TCGA  Populations '

# pop_plot <- ggplot(data = pca, aes_string(paste0('PC',a) , paste0('PC',b), col = paste(cluster) )) + 
#   geom_point(shape = 10, size = 3) +  scale_color_brewer(palette = color, na.value = 'black') +
#   #geom_point(shape = 10, size = 3) +  scale_fill_manual(palette = color_blind, na.value = 'black') +
#   xlab(paste0('PC', a, ' (', signif(pve$p.v.e[a], 3), '%)')) +
#   ylab(paste0('PC', b, ' (', signif(pve$p.v.e[b], 3), '%)')) +
#   labs(color = legend, subtitle = plot_subtitle) + 
#   ggtitle(plot_title) + 
#   theme_update(plot.title = element_text(hjust = 0.5, size = 15)) +
#   theme(axis.title.x = element_text(size = 12, vjust = 1),
#         axis.title.y = element_text(size = 12, vjust = 1),
#         plot.subtitle = element_text(hjust = 0.5),
#         panel.background = element_rect(fill = "white",
#                                         colour = "#F9F9F9",
#                                         size = 0.5, linetype = "solid"))

# ggsave(plot = pop_plot, path = paste0(plink_prefix, '/'),
#        height = 1784, width = 3171, unit = 'px', dpi = 300,
#        filename = paste0(file_pre, '_', cluster, '_PC', a, '_PC', b, '.png'),
#        device = 'png')


# }