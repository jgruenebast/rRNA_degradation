#!/usr/bin/perl -w

use strict;

open(GFF, "<", "/local/projects-t3/SerreDLab-3/jgruenebast/Genomes/Pf_3D7/PlasmoDB-63_Pfalciparum3D7.gff");

my %coding;

while (<GFF>) {
    chomp;

    my $gffinfo = $_;
    my @gene = split /\t/, $gffinfo;

    if ((defined $gene[8]) and ($gene[2] eq "protein_coding_gene")) {
        for (my $bp = $gene[3]; $bp < $gene[4]; $bp++) {
            my $position = "$gene[0]_$gene[6]_$bp";
            $coding{$position} = $gene[8];
        }
    }
}

close GFF;

my $file = "";
for (my $f = 0; $f < 8; $f++) {

	if ($f == 0) {$file = "GamV"}
	if ($f == 1) {$file = "Trophs"}
	if ($f == 2) {$file = "GamVpolyA_new"}
	if ($f == 3) {$file = "TrophspolyA_new"}
	if ($f == 4) {$file = "GamVJHU_new"}
	if ($f == 5) {$file = "Sporozoites"}
	if ($f == 6) {$file = "Rings"}
	if ($f == 7) {$file = "Schizont"}
	print "$file\n";
	
	open (BAM, "samtools view /local/projects-t3/SerreDLab-3/jgruenebast/ONT/NF54/bam-files_NF54/ONT_$file\.bam |");
	open (LEN, ">", "/local/projects-t3/SerreDLab-3/jgruenebast/ONT/rRNA/mRNA_degradation/NF54_$file\_mRNA-all\.csv");
	
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
