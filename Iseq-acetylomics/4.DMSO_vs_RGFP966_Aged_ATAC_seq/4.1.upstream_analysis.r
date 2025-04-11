####质控 比对 去重

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行环境：base

#运行位置：/data/02.private/guanpb/H4K77ac/ATAC

###3例RGFP966处理的 4例DMSO处理的

###原始Fastq存储位置：/home/hulab/K77_rawdata/CUT_tag/ABFC20221500-33/20240627_LH00524_0068_A227VJVLT4

#!/bin/bash

mkdir fastp_report
mkdir cleandata
RAWDATA_FILE="/home/hulab/K77_rawdata/CUT_tag/ABFC20221500-33/20240627_LH00524_0068_A227VJVLT4"

SAMPLE="0uM-1 0uM-2 0uM-3 0uM-4 1uM-1 1uM-2 1uM-4"
for name in $SAMPLE
do
	fileR1="${RAWDATA_FILE}/ATAC-seq_${name}*.R1.fastq.gz" 
	fileR2="${RAWDATA_FILE}/ATAC-seq_${name}*.R2.fastq.gz" 
	fastp -i ${fileR1} -o cleandata/${name}_R1_clean_data.fastq.gz \
	-I ${fileR2} -O cleandata/${name}_R2_clean_data.fastq.gz \
	--thread=5 --html fastp_report/${name}_fastp.html --json fastp_report/${name}_fastp.json
done


SAMPLE="0uM-1 0uM-2 0uM-3 0uM-4 1uM-1 1uM-2 1uM-4"
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

SAMPLE_NAME="0uM-1 0uM-2 0uM-3 0uM-4 1uM-1 1uM-2 1uM-4"

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

####callpeaks

mkdir MACS2_callpeak
cd MACS2_callpeak
bamfile="/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam"

SAMPLE_NAME="0uM-1 0uM-2 0uM-3 0uM-4 1uM-1 1uM-2 1uM-4"
for name in $SAMPLE_NAME
do
	macs2 callpeak \
	-n ${name} \
	-g mm --keep-dup all \
	--treatment ${bamfile}/${name}.sorted.markdup.bam \
	-B --SPMR \
	--call-summits \
	--outdir ${name}_peaks --shift -100 --extsize 200
done

cd /data/02.private/guanpb/H4K77ac/ATAC

SAMPLE_NAME="0uM-1 0uM-2 0uM-3 0uM-4 1uM-1 1uM-2 1uM-4"
for name in $SAMPLE_NAME
do
	bedtools intersect -a /data/02.private/guanpb/H4K77ac/ATAC/MACS2_callpeak/${name}_peaks/${name}_peaks.narrowPeak -b '/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed' -v > ${name}_peaks.rmbl.bed
done


##########合并样本分析
#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行环境：base

#运行位置：/data/02.private/guanpb/H4K77ac/ATAC

#!/bin/bash

mkdir -p mergebam/rmbam

samtools merge mergebam/K77_ATAC_0um.bam \
/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam/0uM-1.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam/0uM-2.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam/0uM-3.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam/0uM-4.sorted.markdup.bam



samtools merge mergebam/K77_ATAC_1um.bam \
/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam/1uM-1.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam/1uM-2.sorted.markdup.bam \
/data/02.private/guanpb/H4K77ac/ATAC/rmdup/rmdupbam/1uM-4.sorted.markdup.bam 


SAMPLE_NAME='K77_ATAC_0um K77_ATAC_1um'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done

cd mergebam/rmbam/
SAMPLE_NAME='K77_ATAC_0um K77_ATAC_1um'

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

cd /data/02.private/guanpb/H4K77ac/ATAC
SAMPLE_NAME='K77_ATAC_0um K77_ATAC_1um'
mkdir bgfile
for name in $SAMPLE_NAME
do
	bedtools genomecov -bg -ibam /data/02.private/guanpb/H4K77ac/ATAC/mergebam/rmbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done

####callpeaks
#!/bin/bash
cd MACS2_callpeak
bamfile="/data/02.private/guanpb/H4K77ac/ATAC/mergebam/rmbam"

SAMPLE_NAME="K77_ATAC_0um K77_ATAC_1um"
for name in $SAMPLE_NAME
do
	macs2 callpeak \
	-n ${name} \
	-g mm --keep-dup all \
	--treatment ${bamfile}/${name}.sorted.markdup.bam \
	-B --SPMR \
	--call-summits \
	--outdir ${name}_peaks --shift -100 --extsize 200
done

SAMPLE_NAME="K77_ATAC_0um K77_ATAC_1um"
for name in $SAMPLE_NAME
do
	bedtools intersect -a /data/02.private/guanpb/H4K77ac/ATAC/MACS2_callpeak/${name}_peaks/${name}_peaks.narrowPeak -b '/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed' -v > rmbl/${name}_peaks.rmbl.bed
done

