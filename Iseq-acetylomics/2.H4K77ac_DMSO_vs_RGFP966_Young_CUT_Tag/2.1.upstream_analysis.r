#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行环境：base

#运行位置：/home/hulab/RGFP966

###3例RGFP966处理的 4例DMSO处理的

###原始Fastq存储位置：/data2/220311_X517_0798_AHHCNGCCX2

#!/bin/bash
mkdir fastp_report
mkdir cleandata
RAWDATA_FILE="/data2/220311_X517_0798_AHHCNGCCX2"
SAMPLE="B1-K77-0 B2-K77-0 B3-K77-0 B4-K77-0 B1-K77-1 B2-K77-1 B3-K77-1 K27-0 K27-1"
for name in $SAMPLE
do
	fileR1="${RAWDATA_FILE}/${name}*.R1.fastq.gz" 
	fileR2="${RAWDATA_FILE}/${name}*.R2.fastq.gz" 
	fastp -i ${fileR1} -o cleandata/${name}_R1_clean_data.fastq.gz \
	-I ${fileR2} -O cleandata/${name}_R2_clean_data.fastq.gz \
	--thread=5 --html fastp_report/${name}_fastp.html --json fastp_report/${name}_fastp.json
done



mkdir -p alignment/sam/bowtie2_summary
mkdir -p alignment/bam
mkdir -p alignment/bed
mkdir -p rmdup/sortbam
mkdir -p rmdup/FixandPos
mkdir -p rmdup/rmdupbam
mkdir -p rmdup/rmdupbam/rmdup_summary



SAMPLE="B1-K77-0 B2-K77-0 B3-K77-0 B4-K77-0 B1-K77-1 B2-K77-1 B3-K77-1 K27-0 K27-1"


for name in $SAMPLE
do	
	bowtie2 --very-sensitive-local --no-mixed --no-discordant -X 1000 -x /data4/cut_TAG/H4K77ac/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
	-1 cleandata/${name}_R1_clean_data.fastq.gz \
	-2 cleandata/${name}_R2_clean_data.fastq.gz -p 8 -S alignment/sam/${name}.sam &> alignment/sam/bowtie2_summary/${name}_bowtie2.txt
	samtools view -o alignment/bam/${name}.bam alignment/sam/${name}.sam
	samtools sort -n alignment/bam/${name}.bam -o rmdup/sortbam/${name}.sorted.bam
	samtools fixmate -m rmdup/sortbam/${name}.sorted.bam rmdup/FixandPos/${name}.fixmate.bam
	samtools sort  rmdup/FixandPos/${name}.fixmate.bam -o rmdup/FixandPos/${name}.positionsort.bam
	samtools markdup -r rmdup/FixandPos/${name}.positionsort.bam rmdup/rmdupbam/${name}.sorted.markdup.bam -s -f rmdup/rmdupbam/rmdup_summary/${name}.rmdup.txt
	
done


#!/bin/bash
cd rmdup/rmdupbam/
SAMPLE="B1-K77-0 B2-K77-0 B3-K77-0 B4-K77-0 B1-K77-1 B2-K77-1 B3-K77-1 K27-0 K27-1"
for name in $SAMPLE
do
	samtools index ${name}.sorted.markdup.bam
done

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！
##运行环境：base
#运行位置：/home/hulab/RGFP966

#!/bin/bash

mkdir -p mergebam/rmbam

samtools merge mergebam/K77-0.bam \
rmdup/rmdupbam/B1-K77-0.sorted.markdup.bam \
rmdup/rmdupbam/B2-K77-0.sorted.markdup.bam \
rmdup/rmdupbam/B3-K77-0.sorted.markdup.bam \
rmdup/rmdupbam/B4-K77-0.sorted.markdup.bam 
 
samtools merge mergebam/K77-1.bam \
rmdup/rmdupbam/B1-K77-1.sorted.markdup.bam \
rmdup/rmdupbam/B2-K77-1.sorted.markdup.bam \
rmdup/rmdupbam/B3-K77-1.sorted.markdup.bam 


SAMPLE_NAME='K77-0 K77-1'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done

cd mergebam/rmbam/
SAMPLE_NAME='K77-0 K77-1'

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

cd /home/hulab/RGFP966/
SAMPLE_NAME='K77-0 K77-1'
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'

mkdir bgfile
for name in $SAMPLE_NAME
do
	bedtools genomecov -bg -ibam /home/hulab/RGFP966/mergebam/rmbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done


SAMPLE_NAME='K77-0 K77-1'
mkdir SEACR_005

SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"

for name in $SAMPLE_NAME
do
	bash $SEACR bgfile/${name}.bedGraph 0.05 non stringent SEACR_005/${name}_005
	bedtools intersect -a SEACR_005/${name}_005.stringent.bed -b ${mmbt} -v >SEACR_005/${name}_005.rmbl.bed
done

#!/bin/bash

mkdir -p mergebam/rmbam

samtools merge mergebam/K77-0_NO_B4.bam \
rmdup/rmdupbam/B1-K77-0.sorted.markdup.bam \
rmdup/rmdupbam/B2-K77-0.sorted.markdup.bam \
rmdup/rmdupbam/B3-K77-0.sorted.markdup.bam  
 


SAMPLE_NAME='K77-0_NO_B4'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done

cd mergebam/rmbam/
SAMPLE_NAME='K77-0_NO_B4'

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

cd /home/hulab/RGFP966/
SAMPLE_NAME='K77-0_NO_B4'
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'

mkdir bgfile
for name in $SAMPLE_NAME
do
	bedtools genomecov -bg -ibam /home/hulab/RGFP966/mergebam/rmbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done


########因为DMSO组有4个样本 其中一个相关性较低 并且B1-B3是一个人做的 B4是另外一个人做的 因此下游分析使用的是没有B4的合并结果

SAMPLE_NAME='K77-0_NO_B4'
mkdir SEACR_005

SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"

for name in $SAMPLE_NAME
do
	bash $SEACR bgfile/${name}.bedGraph 0.05 non stringent SEACR_005/${name}_005
	bedtools intersect -a SEACR_005/${name}_005.stringent.bed -b ${mmbt} -v >SEACR_005/${name}_005.rmbl.bed
done

####由于20250205 revised 版本没有继续使用这个IGV了
######IGV可视化文件位置：/home/hulab/RGFP966/bgfile
##DMSO组：/home/hulab/RGFP966/bgfile/K77-0_NO_B4.bedGraph
##RGFP966组:/home/hulab/RGFP966/bgfile/K77-1.bedGraph
