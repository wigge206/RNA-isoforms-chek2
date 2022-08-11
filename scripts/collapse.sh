## Local paths for files and flair 
Hadley_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/Hadley ## 5 sample fastq.gz files
jessie_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/jessieRun10.fastq.gz
hg38=/STORAGE/Genomes/Gencode/hg38/GRCh38.primary_assembly.genome.fa
gtf=/STORAGE/Genomes/Gencode/hg38/gencode.v41.annotation.gtf
flair=/WORKSPACE/George/Tools/flair/flair.py

mkdir collapse
cd collapse

## Jessie's samples
python $flair collapse -q ../correct/jessieRun10_all_corrected.bed -g $hg38 -r $jessie_fq -t 10  -f $gtf -o jessieRun10_all

## Hadley'samples
python $flair collapse -q ../correct/barcode91_95.bed -g $hg38 -r $hadley_fqPath/* -t 10  -f $gtf -o barcode91_95
