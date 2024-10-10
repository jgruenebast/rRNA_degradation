#!/usr/bin/perl -w

###############################################################################
# This code counts the total number of reads (with a 80% overlap) at:         #
# protein-coding genes, rRNAs, annotated ncRNAs and pseudogenes               #
# Input 1: gff file                                                           #
# Input 2: bam file                                                           #
# Output: Table with total read counts at categories for each file            #   
###############################################################################

use strict;

#open gff file
open(GFF, "<", "/path/to/gff/PlasmoDB-63_Pfalciparum3D7.gff"); #change file path

#define four hash to store the reads for different categories
#4 categories: protein coding genes, rRNAs, pseudogenes, annotated ncRNAs
my %protein_coding_genes;
my %rRNA_genes;
my %pseudogenes;
my %ncRNA_genes;

while (<GFF>) {
    chomp;
    
    #split the table by tab
    my @gene = split /\t/, $_;
    #$gene[0] = Chromosome
		#$gene[2] = protein_coding_gene or rRNA or ncRNA -> check structure of gff file for correct annotations
		#$gene[3] = start position of the gene
		#$gene[4] = end position of the gene
		#$gene[6] = strand
		#$gene[8] = GeneID and description

    #print to make sure that splitting of the gff worked
    #print "$gene[0]\t#$gene[2]\n";

    #store the positions of the protein coding genes, rDNA genes, ncRNAs and pseudogenes from the gff file in a hash
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
#adjust according to numbers of files
#copy the correct file name from folder
#change file names and file path for Dd2 files
my $file = "";
for (my $f = 0; $f < 8; $f++) { 

	if ($f == 0) {$file = "Rings"}
	if ($f == 1) {$file = "Trophozoites"}
	if ($f == 2) {$file = "Schizonts"}
	if ($f == 3) {$file = "GamV"}
	if ($f == 4) {$file = "GamV_2"}
	if ($f == 5) {$file = "Sporozoites"}
	if ($f == 6) {$file = "Trophs_polyA"}
	if ($f == 7) {$file = "GamV_polyA"}
	print "$file\n";
	
	#path to bam files and output files
	open (BAM, "samtools view /path/to/bam/$file\.bam |"); #Change file path to bam files (they need to be in one folder)
	open (OUT, ">", "/path/to/output/$file\.csv"); #Change file path to output files, modify output name if needed eg. NF54_$file\.csv

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

      #column[1] = strand (0 or 16)
      #column [2] = Chromosome
      #column [3] = start position of the aligned read
      #column [5] = Cigar string
      #column [9] = sequence
		
		  #correct the strand annotation for bam to match gff file
    	my $strand;
		  if ($column[1] == 16){
			$strand = "-";
		  } elsif ($column[1] == 0){
			$strand = "+";
		  }
		
		  #Define read length by taking the length of the sequence
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

		  #define inside the loop and set to 0, so that it resets the counter after every read
      #counts the number of bp overlapping with the position stored in the hash for the different categories
		  my $protein_coding_overlap = 0;
		  my $rRNA_overlap = 0;
		  my $pseudogene_overlap = 0;
		  my $ncRNA_overlap = 0;
		
		  #compare the reads from the bam file with the four categories from the gff file (by position and basepair)
    	for (my $i = $column[3]; $i < $end; $i++) {
        	my $pos_bam = join("_", $column[2], $strand, $i); #position needs to be the same structure as from the gff file

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
    
    	#If the read overlaps more than 80% with a categorie, its counted. All other reads are counted in other.
      #Adjust the number if different percentage of overlap is supposed to be counted (e.g. 0.5 for 50%)
    	if ($protein_coding_overlap>(0.8 * $readlength_clipped)) {$read_protein_coding++}
    	elsif ($rRNA_overlap>(0.8 * $readlength_clipped)) {$read_rRNA++}
    	elsif ($pseudogene_overlap>(0.8 * $readlength_clipped)) {$read_pseudogene++}
    	elsif ($ncRNA_overlap>(0.8 * $readlength_clipped)) {$read_ncRNA++}
    	else {$read_other++}
	}
	close BAM;

	# Print the results in output file split by tab (change for "," to get a true .csv file)
	print OUT "$file\tProtein coding genes\t$read_protein_coding\n";
	print OUT "$file\trRNA\t$read_rRNA\n";
	print OUT "$file\tpseudogene\t$read_pseudogene\n";
	print OUT "$file\tncRNA_annotated\t$read_ncRNA\n";
	print OUT "$file\tOther\t$read_other\n";

	close OUT;
}

exit;
