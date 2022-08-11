## Local paths for files and flair 
Hadley_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/Hadley ## 5 sample fastq.gz files
jessie_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/jessieRun10.fastq.gz
hg38=/STORAGE/Genomes/Gencode/hg38/GRCh38.primary_assembly.genome.fa
gtf=/STORAGE/Genomes/Gencode/hg38/gencode.v41.annotation.gtf
flair=/WORKSPACE/George/Tools/flair/flair.py

## Path for analysis
path=/WORKSPACE/George/CHEK2_RNAisoforms/
cd $path

mkdir correct
cd correct

## Can loop through both Jessie's and Hadley's data
align=$path"/align/"*.bed
for bed in $align
do
  out=`basename "$bed"`
  out="${out%.*}"
  python $flair correct -q $bed -g $hg38 -f $gtf -o $out -t 8
done

## Check this - merged hadley samples
cat barcode9[0-5]*corrected.bed > barcode91_95.bed
