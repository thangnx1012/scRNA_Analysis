#1. Download rawdata, QC and Triming
export PATH=/media/hkh/8TB/XUANTHANG/TOOLS/sratoolkit.2.10.9-ubuntu64/bin:$PATH

RawData=/media/hkh/8TB/XUANTHANG/StemcellDiff/hESC_nicotine/RawData
Results=/media/hkh/8TB/XUANTHANG/StemcellDiff/hESC_nicotine/Results
mkdir -p $RawData/FastQC $RawData/sra $RawData/Fastq
mkdir -p $Results/STARsolo


# scRNAseq PAIR-END


prefetch --option-file $RawData/*SRR_Acc_List.txt --output-directory $RawData/sra
fastq-dump --split-files $RawData/sra/SRR*/SRR* --outdir $RawData/Fastq --gzip
fastqc -t 6 -o $RawData/FastQC --noextract -f fastq $RawData/Fastq/SRR*
multiqc $RawData $RawData/FastQC/. 

#2 Alignment (USING UCSC Human hg38)
# RNA-seq (Pair-end, Using STARsolo)


# The 8bp file is the index read (fastq_1), the 26bp is the cellular barcode + UMI read (fastq_2) and the 98bp one is the cDNA (fastq_3)


# Using 10X Genomic V1: https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt
# Using 10X Genomic V2: https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
# Using 10X Genomic V3: https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
      --soloCBwhitelist None # if no whitelist is provided.
      --soloCBlen 16 \ 
      --soloCBstart 1 \ #length cell barcode from 1-16
      --soloUMIlen 10 \ # default V2, length=12 for V3
      --soloUMIstart 17 \ #length of UMI from 17-26 (V2) or 19-30 (V3)
      --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode read (cell+UMI)
      --readFilesCommand zcat \ if gzip fastq file
      --genomeLoad LoadAndKeep
      --soloOutFileNames Results/STARsolo/Solo.out/${i}\ 
      --soloCellFilter  EmptyDrops_CR \ #Filter drop 
      --soloFeatures Gene GeneFull SJ Velocyto \ #Gene: mRNA/ #GeneFull: mRNA+nascent-mRNA/ #SJ: splice junctions/ #Velocyto: calculate Spliced, Unspliced, and Ambiguous counts

for i in $RawData/Fastq/*_1.fastq.gz ; do
echo $i;
base=`basename -s _1.fastq.gz $i`
STAR  --runThreadN 7 \
      --soloType CB_UMI_Simple \
      --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
      --genomeDir /media/hkh/8TB/XUANTHANG/References/Refrence_Human_Genome/Homo_sapiens_UCSC_hg38/Sequence/STAR_Index\
      --soloCBwhitelist $RawData/737K-august-2016_10X_V2.txt \
      --readFilesIn $RawData/Fastq/${base}_3.fastq.gz $RawData/Fastq/${base}_2.fastq.gz \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --soloOutFileNames ${base} genes.tsv barcodes.tsv matrix.mtx \
      --outFileNamePrefix $Results/STARsolo/${base} \
      --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR 
done




























