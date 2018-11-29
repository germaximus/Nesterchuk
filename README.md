# Nesterchuk_2018

**Prerequisites:**  
Cutadapt 1.19  
TopHat 2.1.1  
Bowtie 2.2.7  
featureCounts 1.6  
blast 2.7.1+  
Raw data, index files, and custom annotation can be found in the GEO [NCBI repository]()

### Preparing genome annotation and index files
Mouse genomic sequences and annotation files (GRCm38.p4) were downloaded from the [NCBI repository](ftp://ftp.ncbi.nih.gov/genomes/M_musculus/). To obtain genome assembly, download fasta files of individual chromosomes and concatenate them in the ascending order (omit mitochondrial chromosome). Since all animals were females, Y chromosome should also be discarded.

| files               | MD5 check sum (unzipped)         | Description                                               |
| ------------------- |:--------------------------------:| ----------------------------------------------------------|
| GRCm38.p4.rna.fa    | b22037bd202465ce84cee9f6409331e8 | RNA in fasta format, coding + noncoding                   |
| GRCm38.p4.genome.fa | 49e0e5a638620d90990be88c81030923 | Genome sequence (nuclear genome only)                     |
| GRCm38.p4.gbk       | adc1125bf6b2c3b5a52414fe2fe98ac7 | RNA in gene bank format, coding + noncoding               |
| GRCm38.p4.gff3      | ab982471b0b29ebde3d966ec2424253f | Genome annotation                                         | 

**Customizing genome annotation**  
Extrachromosomal contigs and annotations were omitted. 'Gnomon' (Predicted) records from gff file were also omitted and only 'RefSeq' (manually curated) left. Perl and R scripts are included in the GitHub repository.   
```bash
Discard_extrachromosomal_annotation.pl GRCm38.p4.gff3 >GRCm38.p4.custom.gff
Discard_gnomon_annotation.pl >GRCm38.p4.Refseq.gff	# automatically takes GRCm38.p4.custom.gff as an input
```
Remove non-coding RNA genes, leave only coding genes with their mRNA, transcript, exon, and CDS children. Fix the gff annotation from previous script by matching gene coordinates with the childern coordinates (occured due to removal of Gnomon features).
```bash
Discard_noncoding_annotation.R
```
**Preparing non-redundant transcript sequences**  
Parse GRCm38.p4.gbk end extract the longest transcript for each gene.   
```bash
mRNA_extractor.pl GRCm38.p4.gbk	#generates temp3 file as output
```
Fill up their UTRs to 100 nt based on the genomic coordinates (if they are shorter). Takes temp3 file from previous step as input. Make sure GRCm38.p4.fa genome reference is present in the same folder.
```bash
mRNA_genome_filler.pl	#generates mRNA_100.fasta containing 20575 transcripts
```
To generate a non-redundant subset of transcripts, run blast all vs all, then process with a custom script  
```bash
makeblastdb -in mRNA_100.fasta -dbtype nucl #building an index
blastn -outfmt 6 -evalue 0.001 -db mRNA_100.fasta -query mRNA_100.fasta -out blast_result.txt
BLASTNprocessor.pl blast_result.txt	#generates mRNA_100uniq.fasta containing 17738 transcripts
```
**Building necessary index files**  
```bash
bowtie2-build GRCm38.p4.genome.fa ./Mouse_indices/NCBI_genome #indexing mouse genome for bowtie2 and Tophat
bowtie-build Mouse_rmtRNA.fa ./Mouse_indices/Mouse_rmtRNA
bowtie-build mRNA_100uniq.fa ./Mouse_indices/mRNA_100uniq
tophat -G GRCm38.p4.Refseq.coding.gff --transcriptome-index ./tophat-2.1.1/Mouse_indices/Refseq_coding ./bowtie2-2.2.7/Mouse_indices/NCBI_genome #Indexing mouse transcriptome for TopHat
```
 ### Illumina sequencing reads mapping
 **Liver Transcriptome analysis (mRNA-seq)** 
```bash
cutadapt -j 10 --overlap 5 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o trimmed.fastq input.fastq #adapter trimming
bowtie -p 36 --un filtered.fastq bowtie-1.2.1.1/Mouse_indices/Mouse_rmtRNA trimmed.fastq >/dev/null #filtering out ribosomal, mitochondrial, tRNA and PhiX reads
#mapping for gene expression estimate
tophat -p 50 --library-type fr-firststrand --transcriptome-index ../tophat-2.1.1/Mouse_indices/Refseq_coding --no-novel-juncs -o ./mRNA/ ../bowtie2-2.2.7/Mouse_indices/NCBI_genome filtered.fastq #mapping to a transcriptome and a genome
featureCounts -g gene -s 2 accepted_hits.bam -a ./tophat-2.1.1/Mouse_indices/Refseq_coding.gff -o feature.counts #counting gene expression
#mapping for coverage depth plots
bowtie -p 36 -v 2 -m 1 –-nofw --max redundant.fastq ../bowtie-1.2.1.1/Mouse_indices/mRNA_100uniq filtered.fastq >uniq.bwt
```
**Liver Ribo-seq**  
```bash
cutadapt -j 10 -u 1 -m 23 -a AGATCGGAAGAGCACACGTCT --discard-untrimmed -o trimmed.fastq input.fastq
bowtie -p 36 --un filtered.fastq ./bowtie-1.2.1.1/Mouse_indices/Mouse_rmtRNA trimmed.fastq >/dev/null
#mapping for gene expression estimate
tophat -p 50 --library-type fr-secondstrand --transcriptome-index ./tophat-2.1.1/Mouse_indices/Refseq_coding --no-novel-juncs -o ./output_folder ./bowtie2-2.2.7/Mouse_indices/NCBI_genome filtered.fastq
featureCounts -g gene -s 1 accepted_hits.bam -a ./tophat-2.1.1/Mouse_indices/Refseq_coding.gff -o feature.counts
#mapping for coverage depth plots
bowtie -p 36 -v 2 -m 1 –-norc --max redundant.fastq ../bowtie-1.2.1.1/Mouse_indices/mRNA_100uniq filtered.fastq >uniq.bwt
```
### Metagene Coverage Profiles
Although this information can be obtained from Ribo-seq and mRNA-seq genomic alignment files, it is much easier to re-align raw reads to the *mRNA-100uniq.fasta* file prepared earlier because aligned reads will have transcript coordinates (discontinious) instead of genomic (broken down into exons).
```bash
Coverage.pl uniq.bwt #make sure mRNA_100uniq.fasta is in the same folder with the script or add full path inside the sript. Input file should be in native bowtie-1 format. mRNA_100uniq.fasta should be converted to unix format if run on the linux server.
Coverge.R #process coverage files, plot and explore
```
<details><summary><b>Ribosome occupancy plot </b></summary>
 
Green, red, and blue tracks are patient and two healthy controls corespondingly. Grey track is the mRNA coverage from one of healthy controls. Transcripts are aligned by start codon (left panel) or stop codon (right panel). 100 nt unto UTR are added on both sides.
 

</details>
