#!/usr/bin/perl -w

use strict;

my $output = "";
for (my $c = 0; $c < 20; $c++) {

	if ($c == 0) {$output = "PF3D7_0406200"}
	if ($c == 1) {$output = "PF3D7_1103100"}
	if ($c == 2) {$output = "PF3D7_0708400"}
	if ($c == 3) {$output = "PF3D7_0210100"}
	if ($c == 4) {$output = "PF3D7_1211800"}
	if ($c == 5) {$output = "PF3D7_0818900"}
	if ($c == 6) {$output = "PF3D7_1016900"}
	if ($c == 7) {$output = "PF3D7_0316800"}
	if ($c == 8) {$output = "PF3D7_1331800"}
	if ($c == 9) {$output = "PF3D7_0309600"}
	if ($c == 10) {$output = "PF3D7_1242700"}
	if ($c == 11) {$output = "PF3D7_1309100"}
	if ($c == 12) {$output = "PF3D7_1402500"}
	if ($c == 13) {$output = "PF3D7_0523000"}
	if ($c == 14) {$output = "PF3D7_0827900"}
	if ($c == 15) {$output = "PF3D7_0219400"}
	if ($c == 16) {$output = "PF3D7_0306300"}
	if ($c == 17) {$output = "PF3D7_0517000"}
	if ($c == 18) {$output = "PF3D7_0817500"}
	if ($c == 19) {$output = "PF3D7_1308300"}
	print "$output\n";
	open (GFF, "<", "/path/to/gff/PlasmoDB-63_Pfalciparum3D7.gff"); #change file path

my %coding;

while (<GFF>) {
	chomp;
	
	my $gffinfo = $_;
	my @gene = split /\t/, $gffinfo;
    #$gene[0] = Chromosome
		#$gene[2] = protein_coding_gene or rRNA or ncRNA -> check structure of gff file for correct annotations
		#$gene[3] = start position of the gene
		#$gene[4] = end position of the gene
		#$gene[6] = strand
		#$gene[8] = GeneID and description
  
    #print to make sure that splitting of the gff worked
    #print "$gene[0]\t#$gene[2]\n";
	if((defined $gene[8]) and ($gene[8] =~ "$output")){
		for(my $bp = ($gene[3]); $bp < ($gene[4]); $bp++){	 
			my $position = "$gene[0]_$gene[6]_$bp";
			$coding{$position} = $gene[8];
		}
	}
}
close GFF;

#create a loop for all files
#adjust according to numbers of files
#copy the correct file name from folder
my $file = "";
for (my $f = 0; $f < 11; $f++) {

	if ($f == 0) {$file = "10hpi"}
	if ($f == 1) {$file = "20hpi"}
	if ($f == 2) {$file = "30hpi"}
	if ($f == 3) {$file = "40hpi"}
	if ($f == 4) {$file = "SexRings"}
	if ($f == 5) {$file = "GamI"}
	if ($f == 6) {$file = "GamII"}
	if ($f == 7) {$file = "GamIII"}
	if ($f == 8) {$file = "GamIV"}
	if ($f == 9) {$file = "GamV"}
	if ($f == 10) {$file = "GamVd16"}
	print "$file\n";
	
	open (BAM, "samtools view /path/to/bamfiles/ONT_$file\.bam |"); #change file path
	open (LEN, ">", "/path/to/output/$file\_$output.csv"); #change file path	

my @counts;
my $number_of_coding;
my $read_length;

while (<BAM>) {
	chomp;
	
	#split the table by tab
	my $read = $_;
	my @column = split /\t/, $read; 	
	
  #column [1] = strand (0 or 16)
  #column [2] = Chromosome
  #column [3] = start position of the aligned read
  #column [5] = Cigar string
  #column [9] = sequence
  #print "$column[1]\n";

  #overwrite read strand
	my $strand;
	if ($column[1] == 16){
		$strand = "-";
	} elsif ($column[1] == 0){
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

  #Set counter to 0
	$number_of_coding = 0;
	for(my $i = $column[3]; $i < $end; $i++){
		my $posbam = join("_", $column[2], $strand, $i);
		if (defined $coding{$posbam}){
			$number_of_coding++;
		}
	}

  # change percentage
	my $math = $number_of_coding/$readlength_clipped;
	if($math >= 0.80){
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
}
exit;
