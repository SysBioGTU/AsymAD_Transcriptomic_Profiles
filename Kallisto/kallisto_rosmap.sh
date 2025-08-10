#!/bin/bash
# To index reference transcriptome

kallisto index -i transcriptome.idx Homo_sapiens.GRCh38.cdna.all.fa



 data=/mnt/i/ROSMAP_AD_166/RESULTS_AD
 output=/mnt/c/Users/User/Desktop/Betul_Kallisto/ROSMAP_AD_Rumeysa

 cd $data

 files=($(dir))
 
 for ((a=0; a<(${#files[@]}); a++))
 do
	 samp=${files[$a]}
	 cd $samp
	 echo ${samp} 
	 f1=$(ls *trimmed_p1.fastq*)
	 f2=$(ls *trimmed_p2.fastq*)
	 echo $f1
	 echo $f2
	 Sn=$(echo ${samp}| cut -d'_' -f 3)
     echo ${Sn}


mkdir $output/$samp 
# Kallisto transc exp
kallisto quant -i /mnt/c/Users/User/Desktop/Betul_Kallisto/transcriptome.idx $f1 $f2 -o "$output/$samp" -t 15 --rf-stranded  # please be aware of strand-specificity

cd $data
   
 done

