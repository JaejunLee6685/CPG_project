configfile: "config_1.yaml"
chr = config['chr']

rule all:
    input:
        expand('out1/s4_gzip/TCGA.{chr}.modified.vcf.gz', chr = chr) # with header changed VCF
        expand('out1/s6_add_TCGA_ID/TCGA.modified.splitted.ID.{{chr}}.vcf.gz', chr = chr) # with gnomAD ID attached VCF

rule s001_extract_header:
    input:
        vcfs = f'/storage/scratch01/users/mmoradiellos/tcga_raw_data_from_j.cell.2018.03.039/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.{{chr}}.vcf' # directory of raw vcf splitted with chr
    output:
        vcf_header = temp(f'out1/s1_extract_header/TCGA.{{chr}}.header.vcf.gz') # without header
    resources:
        mem_mb = 12000,
        runtime = 60
    shell:
        'bcftools view -h -O z {input.vcfs} -o {output.vcf_header}'

rule s002_change_header:
    input:
        vcf_header = f'out1/s1_extract_header/TCGA.{{chr}}.header.vcf.gz'
    output:
        vcf_modified_header = f'out1/s2_modified_header/TCGA.{{chr}}.header.modified.vcf'
    resources:
        mem_mb = 12000,
        runtime = 60
    shell:
        "zcat {input.vcf_header} | sed -e 's/ID=AD,Number=R/ID=AD,Number=./g' > {output.vcf_modified_header}" # change header of AD -> to all value to make multiple AD split

rule s003_reheader:
    input:
        vcf_modified_header = f'out1/s2_modified_header/TCGA.{{chr}}.header.modified.vcf',
        vcfs = f'/storage/scratch01/users/mmoradiellos/tcga_raw_data_from_j.cell.2018.03.039/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.{{chr}}.vcf'
    output:
        vcf_modified = temp(f'out1/s3_reheader/TCGA.{{chr}}.modified.vcf')
    resources:
        mem_mb = 12000,
        runtime = 60
    shell:
        "bcftools reheader -h {input.vcf_modified_header} -o {output.vcf_modified} {input.vcfs}" # vcf reheader with modified header -> can make split again

rule s004_gzip:
    input:
        vcf_modified = temp(f'out1/s3_reheader/TCGA.{{chr}}.modified.vcf')
    output:
        vcf_gzip = f'out1/s4_gzip/TCGA.{{chr}}.modified.vcf.gz'
    resources:
        mem_mb = 12000,
        runtime = 60
    shell:
        "gzip -c {input.vcf_modified} > {output.vcf_gzip}" # compress vcf to vcf.gz -> reduce the memory


rule s005_split_multiallelic:
    input:
        tcga = f'out1/s4_gzip/TCGA.{{chr}}.modified.vcf.gz',
        reference = f'/storage/scratch01/users/jaejlee/003.DATA/001.TCGA_VCF/002.Split_VCF/ucsc_hg19_1-22-X-Y.fa'  # human gene refrence fasta hg19 1-22-X-Y / but Y is not included / highly recommend with fa.fai
    output:
        splitted = temp(f'out1/s5_split_multiallelic/TCGA.modified.splitted.{{chr}}.vcf.gz')
    resources:
        mem_mb = 75000,
        runtime = 180
    shell:
        'bcftools norm -f {input.reference} -m - -O z -o {output.splitted} {input.tcga}'  # multiallelic split -> 0/1 , 0/2  -> {ID1 0/1 0/0 | ID2 0/0 0/1} 


rule s006_ADD_ID:
    input:
        splitted = f'out1/s5_split_multiallelic/TCGA.modified.splitted.{{chr}}.vcf.gz'
    output:
        ids = f'out1/s6_add_TCGA_ID/TCGA.modified.splitted.ID.{{chr}}.vcf.gz'
    resources:
        mem_mb = 75000,
        runtime = 180
    shell:
        "bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -o {output.ids} -O z {input.splitted}" # in VCF ID INFO Field -> fill CHR_POS_REF_ALT -> to follow gnomAD way


# snakemake --use-conda --slurm -j 23 -s 1.VCF_header_change.snakefile -n
