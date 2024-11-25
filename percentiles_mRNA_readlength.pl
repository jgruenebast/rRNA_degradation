#!/usr/bin/perl -w

###############################################################################
# This code calculates the 25th, 50th (median), and 75th percentile read      #
# lengths at all protein coding genes genes                                   #
# Input: gff file and bam file                                  #
###############################################################################

use strict;

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

    #Store the positions of all protein-coding genes (check for correct annotations in gff file)
    if((defined $gene[2]) and ($gene[2] eq "protein_coding_gene")){
        for(my $bp = ($gene[3]); $bp < ($gene[4]); $bp++){     
            my $position = "$gene[0]_$gene[6]_$bp";
            $coding{$position} = $gene[8];
        }
    }
}

close GFF;

open (OUT, ">", "path/to/output/proteincoding_percentiles.csv"); #change file path 

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
    
	  open (BAM, "samtools view /path/to/bam/$file\.bam |"); #Change file path to bam files (they need to be in one folder)
    
    # Hash to store read length counts
    my %read_length_counts = ();

    # Compare reads from the BAM file to the hash %coding from the GFF file
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
        
        my $read_length = length($column[9]);
        
        # Define softclipping numbers from cigar string:
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
        # Check if the correct numbers are extracted: 
        # print "$softclip_start\t$softclip_end\t$cigar\t$read_length\t$readlength_clipped\n";
        # sleep (1);
        
        my $end = $column[3] + $readlength_clipped;

        #set counter to 0
        my $number_of_coding = 0;
        for(my $i = $column[3]; $i < $end; $i++){
            my $posbam = join("_", $column[2], $strand, $i);
            if (defined $coding{$posbam}){
                $number_of_coding++;
            }
        }
        
        # Store read length if the read overlaps with the gene of interest
        my $math = $number_of_coding / $readlength_clipped;
        if($math >= 0.80){
            $read_length_counts{$readlength_clipped}++;
        }
    }
    close BAM;
    
    
	# Check if the hash was populated correctly
	#while (my ($length, $count) = each %read_length_counts) {
    #print "Read length: $length, Count: $count\n";
	#}

    # Create an array with the correct counts for each read length
    my @read_lengths = ();
    while (my ($length, $count) = each %read_length_counts) {
        push @read_lengths, ($length) x $count;
    }

	# Print the read lengths array
    #print join(", ", @read_lengths), "\n";

	# Calculate percentiles
	@read_lengths = sort { $a <=> $b } @read_lengths;

  my $n = @read_lengths;

	my $q1 = $read_lengths[int($n * 0.25)];
	my $q3 = $read_lengths[int($n * 0.75)];
	
	my $median;
	if ($n%2) {
	  $median = $read_lengths[int($n * 0.50)];
	}
	else {
	  $median = ($read_lengths[int($n * 0.50)-1] + $read_lengths[int($n * 0.50)]) /2;
	}

    # Output the percentiles
    print OUT "$file\t";
    print OUT "$q1\t$median\t$q3\n";
}

close OUT;
exit;
