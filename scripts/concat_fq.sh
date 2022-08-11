
outPath=/WORKSPACE/George/CHEK2_RNAisoforms/reads/concatanted_reads/
## Path to Jessie's data 
# -- note I have copied her fastq data from HCS to local for speed - this will be rm after concatanted
path=/WORKSPACE/George/CHEK2_RNAisoforms/reads/Jessie/pass/
cat $path/* > $outPath/jessieRun10.fastq.gz

## explicit paths for Hadley's data.
path=/WORKSPACE/George/Hadley_nanoporeData/20220120_1103_X5_FAQ10897_08340024/fastq_pass/
cd $outPath/Hadley
for file in $path/barcode9*
do
  out=`basename $file`
  cat $file/*fastq.gz > $out".fastq.gz"
done