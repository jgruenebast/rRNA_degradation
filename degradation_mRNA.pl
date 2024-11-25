#!/usr/bin/perl -w

use strict;

open(GFF, "<", "/path/to/gff/PlasmoDB-63_Pfalciparum3D7.gff"); #change file path

my %coding;

while (<GFF>) {
    chomp;

    my $gffinfo = $_;
    my @gene = split /\t/, $gffinfo;
 	
  	#$gene[0] = Chr
	#$gene[2] = protein_coding_gene
	#$gene[3] = start
	#$gene[4] = end
	#$gene[6] = strand
	#$gene[8] = GeneID

     #stores all protein-coding genes, check gff file for correct annotation
    if ((defined $gene[8]) and ($gene[2] eq "protein_coding_gene")) {
        for (my $bp = $gene[3]; $bp < $gene[4]; $bp++) {
            my $position = "$gene[0]_$gene[6]_$bp";
            $coding{$position} = $gene[8];
        }
    }
}

close GFF;

#create a loop for all files
#adjust according to numbers of files
#copy the correct file name from folder
#change file names and file path for Dd2 files
my $file = "";
for (my $f = 0; $f < 5; $f++) { 

	if ($f == 0) {$file = "Rings"}
	if ($f == 1) {$file = "Trophozoites"}
	if ($f == 2) {$file = "Schizonts"}
	if ($f == 3) {$file = "GamV"}
	if ($f == 4) {$file = "Sporozoites"}
	print "$file\n";
	
	#path to bam files and output files
	open (BAM, "samtools view /path/to/bam/$file\.bam |"); #Change file path to bam files (they need to be in one folder)
	open (LEN, ">", "/path/to/output/NF54_$file\_mRNA-all\.csv");
	
my @counts = ();
my $number_of_coding;
my $read_length = ();

while (<BAM>) {
    chomp;

    my $read = $_;
    my @column = split /\t/, $read;

    my $strand;
    if ($column[1] == 16) {
        $strand = "-";
    } elsif ($column[1] == 0) {
        $strand = "+";
    }

    $read_length = length($column[9]);

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
    my $readlength_clipped = $read_length - $softclip_start - $softclip_end;
    #Check if the correct numbers are extracted: 
    #print "$softclip_start\t$softclip_end\t$cigar\t$read_length\t$readlength_clipped\n";
    #sleep (1);
	
	my $end = $column[3] + $readlength_clipped;
	$number_of_coding = 0;
    for (my $i = $column[3]; $i < $end; $i++) {
        my $posbam = join("_", $column[2], $strand, $i);
        if (defined $coding{$posbam}) {
            $number_of_coding++;
        }
    }

	#change percentage if needed 	
    my $math = $number_of_coding / $readlength_clipped;
    if ($math >= 0.8) {
        $counts[$readlength_clipped]++;
    }
}

close BAM;

for (my $i = 0; $i < @counts; $i++) {
	if (defined $counts[$i]) {
	print LEN "$i\t$counts[$i]\n";
	}
}
}
exit;
