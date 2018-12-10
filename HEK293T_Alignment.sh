#####################################
##    Set up Files for Analysis    ##
#####################################

#Using UCSC GENCODE v29 annotations which are for hg38 most recent human genome. 
mkdir file.path
STAR --runMode genomeGenerate --runThreadN 20 --genomeDir file.path --genomeFastaFiles file.path --sjdbGTFfile file.path --sjdbOverhang 29

mkdir file.path
STAR --runMode genomeGenerate --runThreadN 20 --genomeDir file.path --genomeFastaFiles file.path --sjdbGTFfile file.path --sjdbOverhang 39

#Downloaded the 45S pre-RNA for homo sapiens from NCBI (NR_146144.1) as fasta
bowtie-build -f file.path output.path

#Prepare plastid metagene
sort -k1,1 -k4,4n UCSC_Gencode_hg38.gtf > UCSC_Gencode_hg38_sort.gtf 
reformat_transcripts --annotation_files UCSC_Gencode_hg38_sort.gtf --sorted --annotation_format GTF2 --output_format BED UCSC_Gencode_hg38.bed
sort -k1,1 -k2,2n UCSC_Gencode_hg38.bed > UCSC_Gencode_hg38_sorted.bed
bedToBigBed UCSC_Gencode_hg38_sorted.bed path.to.sizes UCSC_Gencode_hg38_sorted.bb

#Run the metagene analysis. Will give lots of overlapping regions because it thinks each transcript is a gene. This is fine for our purposes
mkdir Plastid_Metagene_hg38
metagene generate file.path --annotation_files UCSC_Gencode_hg38_sorted.bb --annotation_format BigBed --landmark cds_start
#Freaked out about gene_id but thats fine. basically useless anyways. 

#####################################
##         Download He Data        ##
#####################################
declare -a arr=("SRR5227448" "SRR5227449" "SRR6327777")
for f in "${arr[@]}"
do
    echo "Working on" $f
    fastq-dump.2.9.0 --gzip $f
done

mv SRR5227448.fastq.gz ./HEK_CHX_1_He.fastq.gz
mv SRR5227449.fastq.gz ./HEK_CHX_2_He.fastq.gz
mv SRR6327777.fastq.gz ./HEK_HRT_He.fastq.gz 


#####################################
##           Pre-Processing        ##
#####################################
#Adapter Trimming based on Truseq Ribo Profile kit. 
function riboprof_basic_align {
    echo Received Input Fastq: "$1".fastq.gz
    echo Using $2 processors
    echo Pre-processing input...
    zcat "$1".fastq.gz | cutadapt -a AGATCGGAAGAGCACACGTCT -j 20 --minimum-length 17 - 2>cutadaptlog_$1.log | bowtie -p $2 -v 0 --un "$1"_unalig.fq file.path - >/dev/null 2>bowtielog_$1.log
    echo "done!"
    
    echo Compressing non-rRNA reads using $2 cores...
    pigz "$1"_unalig.fq -p $2
    
    echo Creating directory ./STARout_$1 and aligning to mm10 using $2 cores...
    mkdir ./STARout_$1
    STAR --runMode alignReads --runThreadN $2 --genomeDir file.path --outFileNamePrefix ./STARout_$1/"$1"_ --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 400 --outFilterMismatchNmax 0 --outFilterMatchNmin 15 --outFilterMultimapNmax 1 --readFilesCommand zcat --readFilesIn "$1"_unalig.fq.gz
    echo "done!"
}

declare -a samples=("HEK_CHX_1_He" "HEK_CHX_2_He" "HEK_HRT_He")

for i in ${samples[@]}
do
    riboprof_basic_align $i 10
done

#Clean up Intermediate files 
mkdir ./logs
mv *log ./logs
rm *unalig* #remove intermediate unaligned reads file. 

#####################################
##           Post-Alignment        ##
#####################################
#Index Everything
for_postproc=()
while IFS=  read -r -d $'\0'; do
    for_postproc+=("$REPLY")
done < <(find . -type f -iname "*Aligned.sortedByCoord.out.bam" -print0)

for i in ${for_postproc[@]}; do samtools index $i; done

#Merge CHX replicates
mkdir merged_CHX
samtools merge --output-fmt BAM --threads 20 file.path.merged path.CHX1 path.CHX2
samtools index file.path.merged

#Small Pipeline for only the two samples described 
declare -a for_postproc=("path.CHX.merged" "path.HRT")

#Pipeline.
for i in ${for_postproc[@]}
do
    echo "Found File: $i"
    my_path=$(echo $i | grep -Po ".*(?=_Aligned.sortedByCoord.out.bam)")
    echo "Running psite script. Outfile Name Set as $my_path"
    psite --count_files $i --countfile_format BAM --min_counts 50 --require_upstream --min_length 17 --max_length 35 path.roi.file $my_path
    echo "done!"
    read -p "Please edit the P Offset file located in folder $my_path, then press enter to continue"
    echo Making Wiggle Files. Outfile Name Set as $my_path. Using offset file $my_path"_p_offsets.txt"
    make_wiggle --count_files $i --countfile_format BAM --min_length 17 --fiveprime_variable --offset $my_path"_p_offsets.txt" -o $my_path
    echo "done!"
done

#Make Bigwigs
for_bw=()
while IFS=  read -r -d $'\0'; do
    for_bw+=("$REPLY")
done < <(find . -type f -iname "*wig" -print0)

echo ${for_bw[@]} #really want to be sure. Good thing you are careful about deleting wigs!

for i in ${for_bw[@]}
do
    wigToBigWig "$i" path.to.chrom.sizes "$i.bw" 
done

#Cleaned up 
for i in ${for_bw[@]}; do rm $i; done

#####################################
##      HEK293T RNA-seq Data       ##
#####################################
#Download the single RNA-seq dataset
fastq-dump.2.9.0 --gzip "SRR1630838"
mv SRR1630838.fastq.gz ./HEK293T_RNAseq.fastq.gz

#Run fastqc b/c they don't describe the RNA-seq at all in the paper or in GEO
mkdir fastqc_rnaseq
fastqc -o ./fastqc_rnaseq/ -t 23 HEK293T_RNAseq.fastq.gz

#No obvious adapter contamination. Just align everything with STAR using your RNA-seq settings and go for it. 
mkdir STARout_HEK_RNASeq
STAR --runMode alignReads --runThreadN 23 --genomeDir path.to.dir --outFileNamePrefix path.to.out --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesCommand zcat --outFilterMultimapNmax 10 --readFilesIn HEK293T_RNAseq.fastq.gz

#QC 
samtools index path.to.bam
make_wiggle --count_files path.to.bam  --countfile_format BAM --min_length 17 --fiveprime -o path.to.out
cd ./STARout_HEK_RNASeq
for i in *.wig; do wigToBigWig "$i" path.to.chrom.sizes $i.bw; done;
rm *wig #clean up