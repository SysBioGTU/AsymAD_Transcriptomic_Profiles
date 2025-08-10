#!/bin/bash
# To index reference transcriptome

kallisto index -i transcriptome.idx Homo_sapiens.GRCh38.cdna.all.fa



 data=/mnt/h/MSBB_FP/MSBB_FP_261
 output=/mnt/c/Users/User/Desktop/Betul_Kallisto/Kallisto_FP_Rumeysa

 cd $data

 files=($(dir))
 
 for ((a=0; a<(${#files[@]}); a++))
 do
	 samp=${files[$a]}
	 cd $samp
	 echo ${samp} 
	 f1=$(ls *trimmed_p1.fastq*)
	 echo $f1
	 Sn=$(echo ${samp}| cut -d'_' -f 3)
     echo ${Sn}


mkdir $output/$samp 
# Kallisto transc exp
kallisto quant -i /mnt/c/Users/User/Desktop/Betul_Kallisto/transcriptome.idx -o "$output/$samp" -t 15 --single -l 150 -s 20 $f1

cd $data
   
 done
