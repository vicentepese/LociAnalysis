#!/bin/bash

chrpre=`jq '.chrArray.getSNP' options.json` 
if [[ -z "$chrpre" ]] ||  [[ $chrpre == *"all"* ]]
then 
    echo "No chromosomes specified. SNP matches will be performed in all the chromosomes"
    cd submits
    sbatch --array="1-22" submit_getTransSNP.sh
else
    chr=()
    for i in $chrpre ; do
        if [[ "$i" != *"["* ]] && [[ "$i" != *"]"* ]]
            then
            chr+=$i
        fi
    done
    echo "The preprocessing will be performed for chromosomes ${chr}"
    cd submits
    sbatch --array="$chr" submit_getTransSNP.sh 
fi