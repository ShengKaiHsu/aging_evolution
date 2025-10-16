#!/bin/bash

##illumina raw sequencing read processing
##SKHsu_2018.10.11

cd /Volumes/Temp2/shengkai/age_CGE/formal_experiment

mkdir ageing_formal_output/
mkdir ageing_formal_output/bam_files
mkdir ageing_formal_output/coverage
mkdir ageing_formal_output/fastq_counts
mkdir ageing_formal_output/final_bam
mkdir ageing_formal_output/Rsubread
mkdir ageing_formal_output/trimmed_bam
mkdir ageing_formal_output/trimmed_fastq

date1=`date`

# ####genome indexing (only needed for the first run)#### 
# echo -e "\n\n\n#######################################"
# echo "step1 genome_indexing"
# echo -e "#######################################"
# #gmap_build -d genome_index -D /Volumes/Temp1/shengkai/reference_Dsim/reference_gsnap/ /Volumes/Temp1/shengkai/reference_Dsim/JMCE01.1.chrom_or_gi_bacteria.fa


# ####mapping####
# echo -e "\n\n\n#######################################"
# echo "step2 read_mapping"
# echo -e "#######################################"

# for i in 518 519 520 521
# 	do
# 	echo -e "\n\n#######################################"
# 	echo "mapping ${i}"
# 	echo -e "#######################################"
# 	#trimming
# #	readtools TrimReads -mottQual 20 -I Pool_${i}a_1.fq.gz -I2 Pool_${i}a_2.fq.gz -O /Volumes/Temp2/shengkai/age_CGE/formal_experiment/ageing_formal_output/trimmed_bam/Pool_${i}_trimmed.bam

# #	readtools ReadsToFastq -I /Volumes/Temp2/shengkai/age_CGE/formal_experiment/ageing_formal_output/trimmed_bam/Pool_${i}_trimmed.bam -O /Volumes/Temp2/shengkai/age_CGE/formal_experiment/ageing_formal_output/trimmed_fastq/Pool_${i}_trimmed

# 	#mapping
# 	gsnap --gunzip -d genome_index -D /Volumes/Temp1/shengkai/reference_Dsim/reference_gsnap/ /Volumes/Temp2/shengkai/age_CGE/formal_experiment/ageing_formal_output/trimmed_fastq/Pool_${i}_trimmed_1.fq.gz /Volumes/Temp2/shengkai/age_CGE/formal_experiment/ageing_formal_output/trimmed_fastq/Pool_${i}_trimmed_2.fq.gz -A sam -k 15 -N 1 -t 12 -m 0.08 -o ageing_formal_output/bam_files/Pool_${i}.mapped.sam 1>> ageing_formal_output/mapping_out.txt 2>> ageing_formal_output/mapping_out.txt

# 	#extracting proper reads
# 	samtools view -h -F 4 -q 20 -b ageing_formal_output/bam_files/Pool_${i}.mapped.sam > ageing_formal_output/bam_files/Pool_${i}.mapped.bam 
# 	sambamba sort ageing_formal_output/bam_files/Pool_${i}.mapped.bam -m 8G -t 8 --tmpdir ./ -o ageing_formal_output/bam_files/Pool_${i}.mapped_sorted.bam
# 	rm ageing_formal_output/bam_files/Pool_${i}.mapped.sam
# done

# echo -e "\n\n\n#######################################"
# echo "Mapping finished"
# echo -e "#######################################"

####demultiplexing####
echo -e "\n\n\n#######################################"
echo "step3 demultiplexing"
echo -e "#######################################"
for i in 519 520 521
	do
	echo -e "\n\n#######################################"
	echo "debarcoding ${i}"
	echo -e "#######################################"
	java -jar -Xmx8G -Dreadtools.barcode_index_delimiter='+' /usr/local/Cellar/readtools/1.5.2/libexec/ReadTools.jar AssignReadGroupByBarcode -bc /Volumes/Temp1/shengkai/barcode_info/barcode_Pool_${i}.txt -I ageing_formal_output/bam_files/Pool_${i}.mapped_sorted.bam  -O ageing_formal_output/final_bam/Pool_${i} --barcodeInReadName TRUE --splitReadGroup true --maximumMismatches 3 --maximumN 2
done

echo -e "\n\n\n#######################################"
echo "step4 Count the reads using Rsubread"
echo -e "#######################################"

ls ageing_formal_output/final_bam/*.bam > age_all_bam.txt
current_dir=`pwd`
echo "${current_dir}/ageing_all_bam.txt"
Rscript /Volumes/Temp1/shengkai/monster/count_Rsubread_Qiagen.R --gtf /Volumes/Temp1/shengkai/reference_Dsim/Genome_annotation/JMCE01.2-exon-only.gtf --is_Paired TRUE --bam age_all_bam.txt --output ageing_formal_output/Rsubread 1>> ageing_formal_output/Rout.txt   2>> ageing_formal_output/Rout.txt

echo -e "\n\n\n#######################################"
echo "step5 Inspect the coverage distribution - magenta genes only"
echo -e "#######################################\n\n\n"
for files in `cat age_all_bam.txt | cut -f 1 -d .`
	do
	echo $files
	mv ${files}.bai ${files}.bam.bai
done

export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=/usr/local/lib/python2.7/dist-packages:$PYTHONPATH
python /usr/local/bin/geneBody_coverage.py -r /Volumes/Temp1/shengkai/reference_Dsim/Genome_annotation/JMCE01.2_long_genes_20p.bed -i age_all_bam.txt -o ./age_Magenta 1> ./ageing_formal_output/ageing_coverage_Magenta_out.txt 2> ./ageing_formal_output/ageing_coverage_Magenta_out.txt &


echo "Start Time"
echo $date1
echo "End Time"
date
