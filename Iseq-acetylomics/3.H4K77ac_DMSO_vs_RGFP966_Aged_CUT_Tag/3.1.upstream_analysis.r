####质控 比对 去重

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行环境：base

#运行位置：/data/02.private/guanpb/H4K77ac/CUTTag

###3例RGFP966处理的 3例DMSO处理的

###原始Fastq存储位置：/home/hulab/K77_rawdata/CUT_tag/20240628_LH00308_0170_A227VVLLT4

#!/bin/bash

mkdir fastp_report
mkdir cleandata
RAWDATA_FILE="/home/hulab/K77_rawdata/CUT_tag/20240628_LH00308_0170_A227VVLLT4"
SAMPLE="0uM-501-X-701 0uM-501-X-702 0uM-501-X-703 1uM-502-X-701 1uM-502-X-702 1uM-502-X-703"
for name in $SAMPLE
do
	fileR1="${RAWDATA_FILE}/H4K77ac_${name}*.R1.fastq.gz" 
	fileR2="${RAWDATA_FILE}/H4K77ac${name}*.R2.fastq.gz" 
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

SAMPLE="0uM-501-X-701 0uM-501-X-702 0uM-501-X-703 1uM-502-X-701 1uM-502-X-702 1uM-502-X-703"
for name in $SAMPLE
do	
	bowtie2 --very-sensitive-local --no-mixed --no-discordant -X 1000 -x /data4/cut_TAG/H4K77ac/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
	-1 cleandata/${name}_R1_clean_data.fastq.gz \
	-2 cleandata/${name}_R2_clean_data.fastq.gz -p 8 -S alignment/sam/${name}.sam &> alignment/sam/bowtie2_summary/${name}_bowtie2.txt
	samtools view -o alignment/bam/${name}.bam alignment/sam/${name}.sam
	rm alignment/sam/${name}.sam
	samtools sort -n alignment/bam/${name}.bam -o rmdup/sortbam/${name}.sorted.bam
	samtools fixmate -m rmdup/sortbam/${name}.sorted.bam rmdup/FixandPos/${name}.fixmate.bam
	samtools sort  rmdup/FixandPos/${name}.fixmate.bam -o rmdup/FixandPos/${name}.positionsort.bam
	samtools markdup -r rmdup/FixandPos/${name}.positionsort.bam rmdup/rmdupbam/${name}.sorted.markdup.bam -s -f rmdup/rmdupbam/rmdup_summary/${name}.rmdup.txt	
done

#########callpeaks

#!/bin/bash

SAMPLE_NAME="0uM-501-X-701 0uM-501-X-702 0uM-501-X-703 1uM-502-X-701 1uM-502-X-702 1uM-502-X-703"
mkdir bgfile
for name in $SAMPLE_NAME
do
        bedtools genomecov -bg -ibam /data/02.private/guanpb/H4K77ac/CUTTag/rmdup/rmdupbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done


mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"
mkdir SEACR_005
for name in $SAMPLE_NAME
do
        bash $SEACR /data/02.private/guanpb/H4K77ac/CUTTag//bgfile/${name}.bedGraph 0.05 non stringent SEACR_005/${name}_005
        bedtools intersect -a SEACR_005/${name}_005.stringent.bed -b ${mmbt} -v >SEACR_005/${name}_005.rmbl.bed
done


#######涉及文章Fig.5k IGV输入文件的获取#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行环境：base

#运行位置：/data/02.private/guanpb/H4K77ac/CUTTag

#!/bin/bash

mkdir -p mergebam/rmbam

samtools merge mergebam/K77_Ct_0um.bam \
/data/02.private/guanpb/H4K77ac/CUTTag/rmdup/rmdupbam/0uM-501-X-701.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/CUTTag/rmdup/rmdupbam/0uM-501-X-702.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/CUTTag/rmdup/rmdupbam/0uM-501-X-703.sorted.markdup.bam

samtools merge mergebam/K77_Ct_1um.bam \
/data/02.private/guanpb/H4K77ac/CUTTag/rmdup/rmdupbam/1uM-502-X-701.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/CUTTag/rmdup/rmdupbam/1uM-502-X-702.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/CUTTag/rmdup/rmdupbam/1uM-502-X-703.sorted.markdup.bam


SAMPLE_NAME='K77_Ct_0um K77_Ct_1um'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done

cd mergebam/rmbam/
SAMPLE_NAME='K77_Ct_0um K77_Ct_1um'

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

cd /data/02.private/guanpb/H4K77ac/CUTTag
SAMPLE_NAME='K77_Ct_0um K77_Ct_1um'
mkdir bgfile
for name in $SAMPLE_NAME
do
	bedtools genomecov -bg -ibam /data/02.private/guanpb/H4K77ac/CUTTag/mergebam/rmbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done
