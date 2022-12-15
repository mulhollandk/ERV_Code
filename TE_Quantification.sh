#!/bin/bash
module load samtools
module load STAR

# - - - Quantification of TE

#trim reads for quality
for prefix in $(ls *.fq | sed -r 's/.fq//' | uniq)
do
~/TrimGalore-0.6.5/trim_galore --path_to_cutadapt ~/.local/bin/cutadapt ${prefix}.fq &
done

#align reads to mouse genome (GRCm38) using STAR
for prefix in $(ls *.fq | sed -r 's/_[12]_trimmed.fq//' | uniq)
do
STAR \
--runThreadN 1 \
--outFileNamePrefix ${prefix} \
--outFilterMultimapNmax 4 \
--outFilterMismatchNmax 2 \
--genomeDir mm9 \
--readFilesIn ${prefix}_1_trimmed.fq ${prefix}_2_trimmed.fq \
--outReadsUnmapped Fastx \
--outSAMattributes NH HI AS nM NM MD jM jI XS
done

#sort sam files
for prefix in $(ls *.sam | sed -r 's/.sam//' | uniq)
do
samtools sort -n ${prefix}.sam -o ${prefix}.sorted.sam
done

#align reads to mouse genes (GRCm38; Ensembl 82) using HTSeq
for prefix in $(ls *.sorted.sam | sed -r 's/.sorted.sam//' | uniq)
do
htseq-count -m intersection-nonempty -o ${prefix}_htseq.sam ${prefix}.sam reference/Mus_musculus.GRCm38.101.gtf > ${prefix}_htseq_counts.txt
do

#prepare sam files for repeatmasker alignment
for prefix in $(ls *_htseq.sam | sed -r 's/_htseq.sam//' | uniq)
do
samtools view -b -T /mm9/*.fa ${prefix}_htseq.sam > ${prefix}_htseq.bam
samtools view ${prefix}_htseq.bam | grep 'no_feature' > ${prefix}_nofeature_htseq.bam
samtools view -b -T mm9.fa ${prefix}_nofeature_htseq.bam > ${prefix}_nofeature_htseq.bam
bedtools bamtofastq -i ${prefix}_nofeature_htseq.bam -fq ${prefix}_nofeature_htseq.fastq
done

for prefix in $(ls *.fastq | sed -r 's/.fastq//' | uniq)
do
STAR \
--runThreadN 1 \
--outFileNamePrefix ${prefix} \
--outFilterMultimapNmax 4 \
--outFilterMismatchNmax 2 \
--genomeDir repeatmasker_reference \
--readFilesIn ${prefix}.fastq \
--outReadsUnmapped Fastx \
--outSAMattributes NH HI AS nM NM MD jM jI XS
done
