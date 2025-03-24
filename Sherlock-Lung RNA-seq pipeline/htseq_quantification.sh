#!/bin/bash
module load samtools/1.19
module load htseq/2.0.4
while getopts ":s:i:o:m:" opt;do
  case $opt in
    s ) sample=$OPTARG;;
    i ) indir=$OPTARG;;
    o ) outdir=$OPTARG;;
    m ) mem=$OPTARG;;
    \? ) echo 'usage: htseq_gencode.sh [-s sample][-i input directory][-o output directory][-m java memory]'
    exit 1;;
  esac
done
shift $(($OPTIND-1))

if [[ -z $mem ]];then 
	mem="60g"
fi 
echo -e "sample:$sample\nindir:$indir\noutdir:$outdir\nmem:$mem\n"

workdir=$outdir/$sample
if [[ ! -d $workdir ]];then
    mkdir -p $workdir
fi

if [[ $sample == TCGA-* ]]||[[ $sample == Y* ]];then
  $stranded="no"
else
  $stranded="yes"
fi


GTF="/fdb/GENCODE/Gencode_human/release_35/gencode.v35.annotation.gtf"

cd $workdir
samtools view -b -F 4 $indir/sorted.bam >$workdir/mapped.bam
samtools sort -n -O BAM $workdir/mapped.bam -o $workdir/mapped.name.bam
htseq-count  -m intersection-nonempty \
-i gene_id \
-r name \
-s $stranded \
$workdir/mapped.name.bam $GTF >$workdir/${sample}_htseq.txt

rm $workdir/mapped.name.bam
rm $workdir/mapped.bam

