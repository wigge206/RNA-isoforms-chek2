## Local paths for files and flair 
Hadley_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/Hadley ## 5 sample fastq.gz files
jessie_fq=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/jessieRun10.fastq.gz
hg38=/STORAGE/Genomes/Gencode/hg38/GRCh38.primary_assembly.genome.fa
gtf=/STORAGE/Genomes/Gencode/hg38/gencode.v41.annotation.gtf
flair=/WORKSPACE/George/Tools/flair/flair.py

## Path for analysis
path=/WORKSPACE/George/CHEK2_RNAisoforms/
cd $path

mkdir quantify
cd quantify
manifest_hadley=$path"sample_manifest_Hadley.tsv"
manifest_jessie=$path"sample_manifest_Jessie.tsv"

python $flair quantify -r $manifest_jessie -i ../collapse/jessieRun10_all.isoforms.fa -o jessieRun10_all_countMatrix.tsv
python $flair quantify -r $manifest_hadley -i ../collapse/barcode91_95.isoforms.fa -o hadley_barcode91_95_countMatrix.tsv