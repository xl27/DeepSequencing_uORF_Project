##############################
##      Prepare  STAR       ##
##############################

#First need to make mm9 STAR directory for alignment. In theory could use ours but want to use the OLD UCSC for mm9 annotations
STAR --runMode genomeGenerate --runThreadN 20 --genomeDir path.to.genome --genomeFastaFiles path.to.fasta --sjdbGTFfile path.to.gtf --sjdbOverhang 29

##############################
##      Download Data       ##
##############################

#Run the downloads 
declare -a arr=("SRR315613" "SRR315614" "SRR315615")
for f in "${arr[@]}"
do
    echo "Working on" $f
    fastq-dump.2.9.0 --gzip $f
done

cat *.fastq.gz > ribprof_ES_HRT90s.fastq.gz

zcat ribprof_ES_HRT90s.fastq.gz | cutadapt -a CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTGAA -j 20 - > /dev/null

#Do the rest programatically. 
for x in {315604..315611}
do 
    z="SRR$x"
    echo "Working on" $z
    fastq-dump.2.9.0 --gzip $z
done

HRT_120=()
for i in {315604..315606}; do HRT_120+=("SRR$i.fastq.gz"); done;
echo "${HRT_120[@]}"    #works!
cat "${HRT_120[@]}" > ribprof_ES_HRT120s.fastq.gz

HRT_150=()
for i in {315607..315609}; do HRT_150+=("SRR$i.fastq.gz"); done;
echo "${HRT_150[@]}"    #works!
cat "${HRT_150[@]}" > ribprof_ES_HRT150s.fastq.gz

HRT_180=()
for i in {315610..315611}; do HRT_180+=("SRR$i.fastq.gz"); done;
echo "${HRT_180[@]}"    #works!
cat "${HRT_180[@]}" > ribprof_ES_HRT180s.fastq.gz

#2nd and third replicate of CHX data
for x in {315624..315627}
do 
    z="SRR$x"
    echo "Working on" $z
    fastq-dump.2.9.0 --gzip $z
done

#Concatenate 2nd CHX replicate (GSM765301)
CHX_2=()
for i in {315624..315626}; do CHX_2+=("SRR$i.fastq.gz"); done;
echo "${CHX_2[@]}"    #works!
cat "${CHX_2[@]}" > ribprof_ES_CHX_2.fastq.gz

#Rename the 3rd CHX replicate (GSM765302)
mv SRR315627.fastq.gz ./ribprof_ES_CHX_3.fastq.gz

############################################################
##      Example Adapter Trimming and rRNA Depletion       ##
############################################################

#Adapter Trimming and rRNA depletion
zcat ribprof_ES_HRT90s.fastq.gz | cutadapt -a CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTGAA -j 10 --minimum-length 17 - 2>cutadaptlog_90s.log | bowtie -p 20 -v 0 --un ribprof_ES_HRT90s_unalig.fq path.to.rRNA - >/dev/null 2>bowtielog_90s.log
gzip ribprof_ES_HRT90s_unalig.fq 

#did this for everything else in the same way

############################################################
##           Example Alignment: Ribo Prof                 ##
############################################################

mkdir ./STARout_CHX_mm9
STAR --runMode alignReads --runThreadN 20 --genomeDir path.to.dir --outFileNamePrefix .output.path --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 400 --outFilterMismatchNmax 0 --outFilterMatchNmin 15 --outFilterMultimapNmax 1 --readFilesCommand zcat --readFilesIn path.to.reads

#Post Analysis
cd ./STARout_CHX_mm9
samtools index ribprof_ES_CHX_mm9Aligned.sortedByCoord.out.bam
psite --count_files ribprof_ES_CHX_mm9Aligned.sortedByCoord.out.bam --countfile_format BAM --min_counts 50 --require_upstream --min_length 17 --max_length 35 path.to.psite Psite_CHX
#edited 29-mer to 14 based on png. otherwise fine.
make_wiggle --count_files ribprof_ES_CHX_mm9Aligned.sortedByCoord.out.bam --countfile_format BAM --min_length 17 --fiveprime_variable --offset Psite_CHX_p_offsets.txt -o riboprof_ES_CHX

#Convert to Bigwigs and remove wigs to save space!
#Make Bigwigs
for i in *.wig
do
    wigToBigWig "$i" mm9.chrom.sizes $i.bw 
done

#Clean Up
rm *.wig

############################################################
##                  RNA-seq Alignment                     ##
############################################################

#Trying RNA-Seq Alignment
#Ran FastQC: reads are 40-43bp with no obvious adapter sequences, overrepresented sequences. 
#Reads are poly-adenylated as per GEO accession- though probably avg fragment >40 so that's why no obvious adapters. 
cutadapt -a "A{100}" -j 20 --minimum-length 17 mRNA_ES.fastq.gz > mRNA_ES_trimmed.fastq
pigz mRNA_ES_trimmed.fastq -p 23

#Generate a new STAR index using a 50bp overhang
STAR --runMode genomeGenerate --runThreadN 23 --genomeDir path.to.dir --genomeFastaFiles path.to.fa --sjdbGTFfile path.to.gtf --sjdbOverhang 49

#Align the RNA-Seq Data
STAR --runMode alignReads --runThreadN 23 --genomeDir path.to.dir --outFileNamePrefix path.to.out --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesCommand zcat --outFilterMultimapNmax 10 --readFilesIn path.to.reads

samtools index output.bam
make_wiggle --count_files output.bam --countfile_format BAM --min_length 17 --fiveprime -o ES_mRNA
for i in *.wig; do wigToBigWig "$i" path.to.chrom.sizes $i.bw; done;
rm *wig #clean up
