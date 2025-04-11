###########涉及Fig.5i的motif分析

##########基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 

##运行位置：/data4/cut_TAG/H4K77ac/use/old/bed/motif

####运行环境: base
#!/bin/bash

SAMPLE="Y O common"
mkdir motif
for name in $SAMPLE
do
	mkdir motif/${name}
	/data1/00.software/00.common/homer/bin/findMotifsGenome.pl /data4/cut_TAG/H4K77ac/use/bed/${name}.great.bed mm10 motif/${name} -size given
done

####输出结果文件位置：
#only young:/data4/cut_TAG/H4K77ac/use/old/bed/motif/Y
#only old :/data4/cut_TAG/H4K77ac/use/old/bed/motif/O
#common:/data4/cut_TAG/H4K77ac/use/old/bed/motif/common
