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
| GRCm38.p4.genome.fa | d74346dac686cc024c087f6f2a2fc3cf | Genome sequence (nuclear genome only)                     |
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
Fill up their UTRs to 100 nt based on the genomic coordinates (if they are shorter). Takes temp3 file from previous step as input. Make sure GRCh38.p12.fna genome reference is present in the same folder.
```bash
mRNA_genome_filler.pl	#generates mRNA_100.fasta containing 20575 transcripts
```
To generate a non-redundant subset of transcripts, run blast all vs all, then process with a custom script  
```bash
makeblastdb -in mRNA_100.fasta -dbtype nucl #building an index
blastn -outfmt 6 -evalue 0.001 -db mRNA_100.fasta -query mRNA_100.fasta -out blast_result.txt
BLASTNprocessor.pl blast_result.txt	#generates mRNA_100uniq.fasta containing 17729 transcripts
```
