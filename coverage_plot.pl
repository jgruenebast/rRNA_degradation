#!/usr/bin/perl -w

#################################################################################
## This Code was written to get a coverage profile of rRNA genes               ##
## It counts how many reads are aligned at the position                        ##
## Output is a table (csv) with the Position and the coverage (number of reads)##
#################################################################################

use strict;

my $output = "";
for (my $c = 0; $c < 4; $c++) {

	if ($c == 0) {$output = "PF3D7_0726000"}
	if ($c == 1) {$output = "PF3D7_0725600"}
	if ($c == 2) {$output = "PF3D7_0532000"}
	if ($c == 3) {$output = "PF3D7_0531600"}
	print "$output\n";
	open(GFF, "<", "/path/to/gff/PlasmoDB-63_Pfalciparum3D7.gff"); #change file path


  my %coding;

  while (<GFF>) {
	  chomp;
	
	  #split the table by tab
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

	  #store the position	
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
  	
    my $number_of_coding;
    my %coverage = ();
    my $read_length = 0;

    # Compare reads from the BAM file to the hash %coding from the GFF file
    while (<BAM>) {
    	chomp;
    	
    	#split the table by tab
    	my $read = $_;
    	my @column = split /\t/, $read; 	
    		#print "$column[1]\n";

      #column[1] = strand (0 or 16)
      #column [2] = Chromosome
      #column [3] = start position of the aligned read
      #column [5] = Cigar string
      #column [9] = sequence
    		
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

      #Define end position of the read based on the clipped read length
    	my $end = $column[3] + $readlength_clipped;

      #Set counter to 0: 
    	$number_of_coding = 0;

    	for(my $i = $column[3]; $i < $end; $i++){
    		my $posbam = join("_", $column[2], $strand, $i);
    		if (defined $coding{$posbam}){
    			$number_of_coding++;
    		}
    	}
    	
    	#Store read length if the read overlaps with the gene of interest
      #Adjust percentage of overlap (here 0%)
    	my $math = $number_of_coding/$readlength_clipped;
    	if ($math >= 0) {
            for (my $i = $column[3]; $i < $end; $i++) {
                my $posbam = join("_", $column[2], $strand, $i);
    				    if (defined $coding{$posbam}){
    				        $coverage{$i}++;
                }
            }
      }
}

# To print the Coverage we have to sort and access the key to %coverage and define it as $pos
for my $pos (sort {$a <=> $b} keys %coverage) {
    print OUT "$pos,$coverage{$pos}\n";
}

close BAM;
close OUT;
}
}

exit;
