## Local paths for files and flair 
Hadley_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/Hadley ## 5 sample fastq.gz files
jessie_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/jessieRun10.fastq.gz
hg38=/STORAGE/Genomes/Gencode/hg38/GRCh38.primary_assembly.genome.fa
gtf=/STORAGE/Genomes/Gencode/hg38/gencode.v41.annotation.gtf
flair=/WORKSPACE/George/Tools/flair/flair.py

## Path for analysis
path=/WORKSPACE/George/CHEK2_RNAisoforms/
cd $path

mkdir align
cd align

## Jessie's sample
python $flair align -g $hg38 -r $jessie_fq -t 8 --quality 30 -v1.3 -o jessieRun10

## Hadley's samples
sampleFastq=$hadley_fqPath/*gz
for file in $sampleFastq
do
  out=`basename "$file"` # basename pulls out the filename from a path
  out="${out%.*}" # removes extension
  out="${out%.*}" # removes extension (Only if gzipped)

python $flair align -g $hg38 -r $file -t 8 -v1.3 --quality 30 -o $out
done