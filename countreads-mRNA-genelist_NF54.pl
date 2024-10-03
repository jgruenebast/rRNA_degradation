#!/usr/bin/perl -w

#####################################################################################################################
#                      												    #
# This code creates a genelist of all genes with read counts for protein-coding genes with reads overlapping >= 80% # 
#                      												    #
#####################################################################################################################

use strict;

#loops through all bam files
my $file = "";
for (my $f = 0; $f < 8; $f++) {

	if ($f == 0) {$file = "GamV_unspliced"}
	if ($f == 1) {$file = "Trophs_unspliced"}
	if ($f == 2) {$file = "GamVpolyA_unspliced"}
	if ($f == 3) {$file = "TrophspolyA_unspliced"}
	if ($f == 4) {$file = "GamVJHU_unspliced"}
	if ($f == 5) {$file = "Sporozoites_unspliced"}
	if ($f == 6) {$file = "Rings_unspliced"}
	if ($f == 7) {$file = "Schizont_unspliced"}
	print "$file\n";
	
	open (BAM, "samtools view /local/projects-t3/SerreDLab-3/jgruenebast/ONT/NF54/bam-files_NF54_3D7-UNspliced/ONT_$file\.bam |");
	open (LEN, ">", "/local/projects-t3/SerreDLab-3/jgruenebast/ONT/rRNA/NF54mRNA-readcount/mRNA_$file\.csv");
	open (GFF, "<", "/local/projects-t3/SerreDLab-3/jgruenebast/Genomes/Pf_3D7/PlasmoDB-63_Pfalciparum3D7.gff");

#create a hash to store the protein-coding genes from gff file (by position and basepair)		
my %coding = ();

#create a hash to count reads per protein-coding gene 
my %gene_counts = ();


while (<GFF>) {
	chomp;
	
	#split the table by tab
	my $gffinfo = $_;
	my @gene = split /\t/, $gffinfo;

		#check if that worked 
		#print "$gene[0]\t$gene[2]\n";
		
		#store each basepair from start to end in a hash
		#$gene[0] = Chr
		#$gene[2] = protein_coding_gene
		#$gene[3] = start
		#$gene[4] = end
		#$gene[6] = strand
		#$gene[8] = GeneID
		
	if((defined $gene[2]) and ($gene[2] eq "protein_coding_gene")){
		for(my $bp = ($gene[3]); $bp < ($gene[4]); $bp++){	 
			my $position = "$gene[0]_$gene[6]_$bp";
			$coding{$position} = $gene[8];
			$gene_counts{$gene[8]} = 0;
			#print "$position\n";
			#sleep(1);
		}
	}
}


my $number_of_coding;
my $gene;

while (<BAM>) {
	chomp;
	
	#split the table by tab
	my $read = $_;
	my @column = split /\t/, $read; 	
		#print "$column[1]\n";
	
 	#correct the strand annotation for bam to match gff file
	my $strand;
	if ($column[1] == 16){
		$strand = "-";
	} elsif ($column[1] == 0){
		$strand = "+";
	}
 
	#define read length	
	my $read_length = length($column[9]);
	
	#Define softclipping numbers from cigar string:
    	my $cigar = $column[5];
    	my $softclip_start = 0; 
    	my $softclip_end = 0;
    	if ($cigar =~ /^(\d+)S/) {
        	$softclip_start = $1;
    	}
    	if ($cigar =~ /(\d+)S$/) {
        $softclip_end = $1;
    	}
     
        #adjust read length based on soft-clipping
    	my $readlength_clipped = $read_length - $softclip_start - $softclip_end;
    	#Check if the correct numbers are extracted: 
    	#print "$softclip_start\t$softclip_end\t$cigar\t$read_length\t$readlength_clipped\n";
    	#sleep (1);

 	#define read end
	my $end = $column[3] + $readlength_clipped;
		#column[1] = strand (0 or 16), column [2] = Chromosome, column [3] =start position
	
 	#set counter to 0 (variable counts how many bp overlap)
	$number_of_coding = 0;
	
	$gene = "";
	for(my $i = $column[3]; $i < $end; $i++){
		my $posbam = join("_", $column[2], $strand, $i);
		if (defined $coding{$posbam}){
			$gene = $coding{$posbam};
			$number_of_coding++;
		}
	}
 
 	#only count reads that overlap more than 80% with a protein coding gene
	my $math = $number_of_coding / $readlength_clipped;
   	if ($math >= 0.8) {
		$gene_counts{$gene}++;
	}

}

close BAM;

#print list with gene names
for my $gene_name (sort {$a cmp $b} keys %gene_counts) {
        print LEN "$gene_name\t$gene_counts{$gene_name}\n";

}
close LEN;
}
close GFF;
exit;
