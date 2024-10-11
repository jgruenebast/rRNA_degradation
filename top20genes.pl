#!/usr/bin/perl -w

use strict;

my $output = "";
for (my $c = 0; $c < 20; $c++) {

	if ($c == 0) {$output = "PF3D7_0406200"}
	if ($c == 1) {$output = "PF3D7_1302100"}
	if ($c == 2) {$output = "PF3D7_0708400"}
	if ($c == 3) {$output = "PF3D7_1246200"}
	if ($c == 4) {$output = "PF3D7_1211800"}
	if ($c == 5) {$output = "PF3D7_0818900"}
	if ($c == 6) {$output = "PF3D7_1016900"}
	if ($c == 7) {$output = "PF3D7_1444800"}
	if ($c == 8) {$output = "PF3D7_0617800"}
	if ($c == 9) {$output = "PF3D7_1105000"}
	if ($c == 10) {$output = "PF3D7_1361800"}
	if ($c == 11) {$output = "PF3D7_1108600"}
	if ($c == 12) {$output = "PF3D7_1003600"}
	if ($c == 13) {$output = "PF3D7_0503400"}
	if ($c == 14) {$output = "PF3D7_0932200"}
	if ($c == 15) {$output = "PF3D7_1248500"}
	if ($c == 16) {$output = "PF3D7_1129100"}
	if ($c == 17) {$output = "PF3D7_1246400"}
	if ($c == 18) {$output = "PF3D7_0818200"}
	if ($c == 19) {$output = "PF3D7_1008700"}
	print "$output\n";
 
	open(GFF, "<", "/path/to/gff/PlasmoDB-63_Pfalciparum3D7.gff"); #change file path

#define a hash to store the positions of the genes
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
#change file names and file path for Dd2 files
my $file = "";
for (my $f = 0; $f < 6; $f++) { 

	if ($f == 0) {$file = "Rings"}
	if ($f == 1) {$file = "Trophozoites"}
	if ($f == 2) {$file = "Schizonts"}
	if ($f == 3) {$file = "GamV"}
	if ($f == 4) {$file = "GamV_2"}
	if ($f == 5) {$file = "Sporozoites"}
	print "$file\n";
	
	open (BAM, "samtools view /path/to/bam/$file\.bam |"); #change file path to bam files (they need to be in one folder)
	open (LEN, ">", "/path/to/output/$file\_$output.csv"); #change file path

  #define an array to store read counts at specific read length
  my @counts;

  my $number_of_coding;
  my $read_length;

  while (<BAM>) {
	  chomp;
	
	  #split the table by tab
	  my $read = $_;
	  my @column = split /\t/, $read; 	
	  #print "$column[1]\n";

    #column [1] = strand (0 or 16)
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
  
      #Define clipped read length
      my $readlength_clipped = $read_length - $softclip_start - $softclip_end;
      #Check if the correct numbers are extracted: 
      #print "$softclip_start\t$softclip_end\t$cigar\t$read_length\t$readlength_clipped\n";
      #sleep (1);
  
      #Define end based on clipped read length
  	  my $end = $column[3] + $readlength_clipped;
  		
  	  $number_of_coding = 0;
  	  for(my $i = $column[3]; $i < $end; $i++){
  		    my $posbam = join("_", $column[2], $strand, $i);
  		    if (defined $coding{$posbam}){
  			      $number_of_coding++;
  		    }
  	  }
     
  	  my $math = $number_of_coding/$readlength_clipped;
     
      #only count reads that overlap more than 80% with a protein coding gene (change if needed)
  	  if($math >= 0.80){
  	     $counts[$readlength_clipped]++;
  	  }
  }

close BAM;

#prints a file with read length in the first column and #reads with that length in the second column
for (my $i = 0; $i < @counts; $i++) {
	if (defined $counts[$i]) {
	print LEN "$i\t$counts[$i]\n";
	}
}
}
}
exit;
