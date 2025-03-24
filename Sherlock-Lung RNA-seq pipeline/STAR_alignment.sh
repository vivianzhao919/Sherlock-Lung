#! /bin/bash
while getopts ":s:t:i:o:g:n:" opt;do
  case $opt in
    s ) sample=$OPTARG;;
    t ) TYPE=$OPTARG;;
    i ) indir=$OPTARG;;
    o ) outdir=$OPTARG;;
    g ) group=$OPTARG;;
    n ) nt=$OPTARG;;
    \? ) echo 'usage: sentieon_rna_star.sh [-s sample][-t raw/trim/plain][-i input directory][-o output directory][-g read group][-n Nthreads]'
    exit 1;;
  esac
done
shift $(($OPTIND-1))
echo -e "sample:$sample\ntype:$TYPE\nindir:$indir\noutdir:$outdir\n"


# Update with the fullpath location of your sample fastq
if [[ $TYPE = "raw" ]];then
	fastq_1=$(ls $indir|grep -E "${sample}(\S)*1.fastq.gz")
	fastq_2=$(ls $indir|grep -E "${sample}(\S)*2.fastq.gz")
	[ -z "$fastq_1" ] && echo "Fastq file of sample $sample is not available"
elif [[ $TYPE = "trim" ]];then
	fastq_1= $(ls $indir|grep -E "${sample}(\S)*1(\S)*trimmed.fastq.gz")
	fastq_2= $(ls $indir|grep -E "${sample}(\S)*2(\S)*trimmed.fastq.gz")
	[ -z "$fastq_1" ] && echo "Fastq file of sample $sample is not available"
elif [[ $TYPE = "plain" ]];then
	echo "plain"
	fastq_1=$(ls $indir|grep -E "${sample}(\S)*1.fastq")
        fastq_2=$(ls $indir|grep -E "${sample}(\S)*2.fastq")
	[ -z "$fastq_1" ] && echo "Fastq file of sample $sample is not available"
else
	echo 'usage: -t raw/trim/plain'
	exit 1
fi

if [[ -z $group ]];then
    group=$sample;
fi
if [[ -z $nt ]];then
    nt=4;
fi
platform="ILLUMINA"


# Other settings
#nt=4 #number of threads to use in computation
workdir=${outdir}/${sample} #Determine where the output files will be stored


if [[ ! -d $workdir ]];then
    mkdir -p $workdir
fi


# Update with the location of the reference data files
fasta="/fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta"
star_fasta="/data/zhaow5/GENOMEDIR/hg38"
GTF="/fdb/GENCODE/Gencode_human/release_35/gencode.v35.annotation.gtf"

# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR=$SENTIEON_INSTALL_DIR
export SENTIEON_LICENSE=$SENTIEON_LICENSE
export SENTIEON_AUTH_DATA=$SENTIEON_LICENSE

star_binary="/usr/local/apps/STAR/2.7.3a/bin/STAR"

module load samtools


# ******************************************
# 0. Setup
# ******************************************

logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir
if [ ! -z "$interval_file" ]; then
  driver_interval_option="--interval $interval_file"
  realign_interval_option="--interval_list $interval_file"
fi

# ******************************************
# 1. Mapping reads with STAR
# ******************************************
if [ -z "$star_fasta" ]; then
  star_fasta="genomeDir"
  mkdir $star_fasta
  $star_binary --runMode genomeGenerate --genomeDir $star_fasta --genomeFastaFiles $fasta --runThreadN $nt
fi
#perform the actual alignment and sorting
$star_binary --twopassMode Basic --genomeDir $star_fasta --runThreadN $nt \
--outSAMtype BAM SortedByCoordinate \
--sjdbGTFfile $GTF \
--twopass1readsN -1 --sjdbOverhang 100 \
--readFilesIn $indir/$fastq_1 $indir/$fastq_2 --readFilesCommand zcat \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattrRGline ID:$sample SM:$sample PL:$platform \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--outSAMattributes NH HI NM MD AS XS nM \
--outSAMunmapped Within \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 8 \
--chimOutJunctionFormat 1 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimMultimapNmax 20 \
--chimNonchimScoreDropMin 10 \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1 \
--alignInsertionFlush Right \
--alignSplicedMateMapLminOverLmate 0 \
--alignSplicedMateMapLmin 30 

mv Aligned.sortedByCoord.out.bam sorted.bam
mv Aligned.toTranscriptome.out.bam Aligned.toTranscriptome.out.ALL.bam
samtools view -b -F 4 Aligned.toTranscriptome.out.ALL.bam >Aligned.toTranscriptome.out.bam

$SENTIEON_INSTALL_DIR/bin/sentieon util index sorted.bam

# ******************************************
# 2. Metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver $driver_interval_option -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

rm -r _STARgenome
rm -r _STARpass1

