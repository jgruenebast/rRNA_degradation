#!/usr/bin/perl -w

###############################################################################
# This code counts the total number of reads (with a 80% overlap) at all:     #
# protein-coding genes, rRNAs, annotated ncRNAs, pseudogenes                  #
###############################################################################

use strict;

#open gff file
open(GFF, "<", "/local/projects-t3/SerreDLab-3/jgruenebast/Genomes/Pf_3D7/PlasmoDB-63_Pfalciparum3D7.gff");

#define hash to store the reads for different categories 
my %protein_coding_genes;
my %rRNA_genes;
my %pseudogenes;
my %ncRNA_genes;

#store the genes from the gff file in a hash (bp by bp) to compare to reads from the bam file
#4 categories: protein coding genes, rRNAs, pseudogenes, annotated ncRNAs
while (<GFF>) {
    chomp;
    
    #split the table by tab
    my @gene = split /\t/, $_;

    if (defined ($gene[2]) and $gene[2] eq 'protein_coding_gene') {
        for (my $bp = $gene[3]; $bp < $gene[4]; $bp++) {
            my $position = "$gene[0]_$gene[6]_$bp";
            $protein_coding_genes{$position} = 1;
        }
    } elsif (defined ($gene[2]) and $gene[2] eq 'rRNA') {
        for (my $bp = $gene[3]; $bp < $gene[4]; $bp++) {
            my $position = "$gene[0]_$gene[6]_$bp";
            $rRNA_genes{$position} = 1;
        }
     } elsif (defined ($gene[2]) and $gene[2] eq 'pseudogene') {
        for (my $bp = $gene[3]; $bp < $gene[4]; $bp++) {
            my $position = "$gene[0]_$gene[6]_$bp";
            $pseudogenes{$position} = 1;
        }
    } elsif (defined ($gene[2]) and $gene[2] eq 'ncRNA_gene') {
        for (my $bp = $gene[3]; $bp < $gene[4]; $bp++) {
            my $position = "$gene[0]_$gene[6]_$bp";
            $ncRNA_genes{$position} = 1;
        }
    }
}
close GFF;

#create a loop for all files
my $file = "";
for (my $f = 0; $f < 14; $f++) {

	if ($f == 0) {$file = "10hpi_DD2_3D7"}
	if ($f == 1) {$file = "20hpi_DD2_3D7_new"}
	if ($f == 2) {$file = "30hpi_2_DD2_3D7"}
	if ($f == 3) {$file = "40hpi_2_DD2_3D7"}
	if ($f == 4) {$file = "SexRings_DD2_3D7_new"}
	if ($f == 5) {$file = "GamI_DD2_3D7_new"}
	if ($f == 6) {$file = "GamII_DD2_3D7_new"}
	if ($f == 7) {$file = "GamIII_DD2_3D7"}
	if ($f == 8) {$file = "GamIV_DD2_3D7_new"}
	if ($f == 9) {$file = "GamV_DD2_3D7"}
	if ($f == 10) {$file = "GamVd16_DD2_3D7"}
	if ($f == 10) {$file = "20hpipolyA_DD2_3D7"}
	if ($f == 12) {$file = "GamVpolyA_DD2_3D7"}
	if ($f == 13) {$file = "GamVd16polyA_DD2_3D7"}
	print "$file\n";
	
	#path to bam files and output files
	open (BAM, "samtools view /local/projects-t3/SerreDLab-3/jgruenebast/ONT/Dd2/bam-files_DD2/3D7Alignment/ONT_$file\.bam |");
	open (OUT, ">", "/local/projects-t3/SerreDLab-3/jgruenebast/ONT/rRNA/countreads80/DD2_$file\.csv");

	#define different categories outside of the loop and set to 0 (this way the total read number is counted)
	my $read_protein_coding = 0;
	my $read_rRNA = 0;
	my $read_pseudogene = 0;
	my $read_ncRNA = 0;
	my $read_other = 0;

	while (<BAM>) {
    	chomp;
    	my $read = $_;
		my @column = split /\t/, $read;
		
		#correct the strand annotation for bam to match gff file
    	my $strand;
		if ($column[1] == 16){
			$strand = "-";
		} elsif ($column[1] == 0){
			$strand = "+";
		}
		
		#Define read length
    	my $read_length = length($column[9]);
    	
		#Define softclipping numbers from cigar string to correct for read length:
    	my $cigar = $column[5];
    	my $softclip_start = 0; 
    	my $softclip_end = 0;
    	if ($cigar =~ /^(\d+)S/) {
        	$softclip_start = $1;
    	}
    	if ($cigar =~ /(\d+)S$/) {
        	$softclip_end = $1;
    	}
    	
    	#Define clipped read length
    	my $readlength_clipped = $read_length - $softclip_start - $softclip_end;
    	
    	#Check if the correct numbers are extracted: 
    	#print "$softclip_start\t$softclip_end\t$cigar\t$read_length\t$readlength_clipped\n";
    	#sleep (1);
    	
		#Define end based on clipped read length
		my $end = $column[3] + $readlength_clipped;

		#define inside the loop so that it resets counter after every read
		my $protein_coding_overlap = 0;
		my $rRNA_overlap = 0;
		my $pseudogene_overlap = 0;
		my $ncRNA_overlap = 0;
		
		#compare the reads from the bam file with the four categories from the gff file (by position and basepair)
    	for (my $i = $column[3]; $i < $end; $i++) {
        	my $pos_bam = join("_", $column[2], $strand, $i);

        	if ($protein_coding_genes{$pos_bam}) {
            	$protein_coding_overlap++;
        	} elsif ($rRNA_genes{$pos_bam}) {
            	$rRNA_overlap++;
        	} elsif ($pseudogenes{$pos_bam}) {
            	$pseudogene_overlap++;
        	} elsif ($ncRNA_genes{$pos_bam}) {
            	$ncRNA_overlap++;
        	} 
    	}
    
    	#If the read overlaps more than 80% to a categorie, its counted. All other reads are counted in other.
    	if ($protein_coding_overlap>(0.8 * $readlength_clipped)) {$read_protein_coding++}
    	elsif ($rRNA_overlap>(0.8 * $readlength_clipped)) {$read_rRNA++}
    	elsif ($pseudogene_overlap>(0.8 * $readlength_clipped)) {$read_pseudogene++}
    	elsif ($ncRNA_overlap>(0.8 * $readlength_clipped)) {$read_ncRNA++}
    	else {$read_other++}
	}
	close BAM;

	# Print the results in output file
	print OUT "$file\tProtein coding genes\t$read_protein_coding\n";
	print OUT "$file\trRNA\t$read_rRNA\n";
	print OUT "$file\tpseudogene\t$read_pseudogene\n";
	print OUT "$file\tncRNA_annotated\t$read_ncRNA\n";
	print OUT "$file\tOther\t$read_other\n";

	close OUT;
}

exit;
