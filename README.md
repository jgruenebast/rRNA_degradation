# rRNA_degradation

Scripts used for analysis of rRNA degradation in _Plasmodium falciparum_ asexual and sexual stages (strain: NF54 and Dd2). 

The data was generated using the SQK-RNA002 kit (Oxford Nanopore Technologies) 

**Step 1:** Reads (fast5 format) were base-called using guppy version 6.4.2 on a GPU 

/usr/local/packages/guppy-6.4.2_gpu/bin/guppy_basecaller -x "cuda:0" --input_path /path/to/fast5 --save_path /path/to/output-directory --config rna_r9.4.1_70bps_hac.cfg --min_qscore 7 --records_per_fastq 10000000 --gpu_runners_per_device 8 --num_callers 1

**Step 2:** Concatenate your fastq files that passed guppy using the cat command

cat *fastq > combined_sample.fastq

**Step 3:** Use minimap for alignment 

minimap2 -ax map-ont -t 2 /path/to/reference_genome.fasta /path/to/guppy/pass/sample.fastq > sample.sam

**Step 4:** Use samtools to sort and convert sam to bam 

samtools view -bhF 2308 sample.sam | samtools sort -o sample.bam

**Step 5:** Index bam file using samtools

samtools index -b sample.bam

**Step 6:** Visualize the data using IGV

**Explanation of different Perl scrips:**

1) Count the number of reads at different categories (mRNA, rRNA, ncRNA, pseudogenes, other): count_totalreads_mRNA-rRNA-ncRNA.pl
2) Count the number or reads at different rRNA types (A, S1, S2): rRNA-types.pl
3) Count the number or full-length reads at different rRNA types (A, S1, S2): rRNA-types_full-length.pl
4) Count consecutive adenosines in the fasta sequence: consecutive_As and convert them into a bedgraph file to look at IGV: covert_As-to-bedgraph.R
5) Count reads to make coverage plots: coverage_plot.pl
6) Count the number of reads of a given read length of reads aligning to rRNA genes to assess degradation: degradation_rRNA
7) Count the number of reads of a given read length of reads aligning to all protein-coding genes: degradation_mRNA.pl
8) Make a list of genes and their read counts: genelist.pl
9) Assess degradation of the 20 most abundant genes in NF54 parasites: top20genes_NF54.pl
10) Assess degradation of the 20 most abundant genes in Dd2 parasites: top20genes_Dd2.pl
11) Calculate the percentiles of read length at rRNAs: percentiles_rRNA_readlength.pl
12) Calculate the percentiles of read length at all mRNAs: percentiles_mRNA_readlength.pl


For SRA upload fast5 files were base-called using an older version of guppy to match SRA criteria and have fast5 files with base-called information. However, all analysis have been based on the above described workflow.  

/usr/local/packages/guppy-4.2.2/bin/guppy_basecaller -x "cuda:6" --input_path /path/to/fast5 --save_path /path/to/output-directory --config rna_r9.4.1_70bps_hac.cfg --gpu_runners_per_device 8 --num_callers 1 --fast5_out
