#1. Download rawdata, QC and Triming
export PATH=/media/hkh/8TB/XUANTHANG/TOOLS/sratoolkit.2.10.9-ubuntu64/bin:$PATH

RawData=/media/hkh/8TB/XUANTHANG/MammaryGland/scRNA_ER/RawData
Results=/media/hkh/8TB/XUANTHANG/MammaryGland/scRNA_ER/Results
mkdir -p $RawData/FastQC $RawData/sra $RawData/Fastq
mkdir -p $Results/STARsolo


# scRNAseq PAIR-END


prefetch --option-file $RawData/SRP180002_SRR_Acc_List.txt --output-directory $RawData/sra
fastq-dump --split-files $RawData/sra/SRR*/SRR* --outdir $RawData/Fastq --gzip
fastqc -t 6 -o $RawData/FastQC --noextract -f fastq $RawData/Fastq/SRR*
multiqc $RawData $RawData/FastQC/. 

#2 Alignment (USING UCSC MOUSE mm10)
# RNA-seq (Pair-end, Using STARsolo)


# The 8bp file is the index read (fastq_3), the 26bp is the cellular barcode + UMI read (fastq_1) and the 101bp one is the cDNA (fastq_2)

# Identify correct cell barcodes in fastq1
# umi_tools whitelist --stdin hgmm_100_R1.fastq.gz \
#                     --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
#                     --set-cell-number=100 \
#                     --log2stderr > whitelist.txt
# Using 10X Genomic V1: https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt
# Using 10X Genomic V2: https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
# Using 10X Genomic V3: https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
      --soloCBwhitelist None # if no whitelist is provided.
      --soloCBlen 16 \
      --soloCBstart 1 \ #length cell barcode from 1-16
      --soloUMIlen 10 \ # default V2, length=12 for V3
      --soloUMIstart 17 \ #length of UMI from 17-26 (V2) or 19-30 (V3)
      --readFilesCommand gzcat \ if gzip fastq file
      --genomeLoad LoadAndKeep
      --soloOutFileNames Results/STARsolo/Solo.out/${i}\
      --soloCellFilter  EmptyDrops_CR \ #Filter drop 
      --soloFeatures Gene GeneFull SJ Velocyto \ #Gene: mRNA/ #GeneFull: mRNA+nascent-mRNA/ #SJ: splice junctions/ #Velocyto: calculate Spliced, Unspliced, and Ambiguous counts

for i in $RawData/sra/*/*.sra; do
echo $i;
base=`basename -s .sra $i`
STAR  --runThreadN 6 \
      --soloType CB_UMI_Simple \
      --genomeDir /media/hkh/8TB/XUANTHANG/References/Reference_Mouse/Mus_musculus_UCSC_mm10/Sequence/STAR_Index \
      --soloCBwhitelist $RawData/737K-august-2016_10X_V2.txt \
      --readFilesIn $RawData/Fastq/${base}_2.fastq $RawData/Fastq/${base}_1.fastq \
      --outSAMtype BAM SortedByCoordinate \
      --soloOutFileNames ${base} genes.tsv barcodes.tsv matrix.mtx \
      --outFileNamePrefix $Results/STARsolo1/${base} \
      --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR 
done




























