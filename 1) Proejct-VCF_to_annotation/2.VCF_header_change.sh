#!/bin/bash

declare -a StringArray=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X')
#declare -a StringArray=('X')
directory_vcf='./out/s6_add_TCGA_ID'
vcfs='TCGA.modified.splitted.ID'
ends='vcf.gz'

for chr in ${StringArray[@]};
do
    echo 'run ${chr}'
    tabix -f -p vcf ${directory_vcf}/${vcfs}.${chr}.${ends} 
done

# this process is to make vcf.gz.tbi -> for the further process with bedtools
# shell script to generate .tbi file -> 