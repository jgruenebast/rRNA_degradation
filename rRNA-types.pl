#!/usr/bin/perl -w
use strict;

###############################################################################
# This code counts the total number of reads (with a 80% overlap) at:         #
# 18S and 28S rRNA types (by Gene_ID)                                         #
# Input 1: gff file                                                           #
# Input 2: bam file                                                           #
# Output: Table with total read counts for each Gene_ID for each file         #   
###############################################################################

my $output = "";
my %coding;

open (OUT, ">", "/path/to/output/rRNA-types.csv"); #change file path

#store all Gene_IDs in a hash
my %gene_ID;
$gene_ID{"PF3D7_0112700"} = 1;
$gene_ID{"PF3D7_0532000"} = 1;
$gene_ID{"PF3D7_0726000"} = 1;
$gene_ID{"PF3D7_1148640"} = 1;
$gene_ID{"PF3D7_1371300"} = 1;
$gene_ID{"PF3D7_0112300"} = 1;
$gene_ID{"PF3D7_0531600"} = 1;
$gene_ID{"PF3D7_0725600"} = 1;
$gene_ID{"PF3D7_1148600"} = 1;
$gene_ID{"PF3D7_1371000"} = 1;

open (GFF, "<", "/path/to/gtf/PlasmoDB-63_Pfalciparum3D7.gff");
while (<GFF>) {
        chomp;
        my $gffinfo = $_;
        my @gene = split /\t/, $gffinfo;
        if((defined $gene[8]) and ($gene[2] eq "ncRNA_gene")) {
        	my @name = split ";", $gene[8];
    		#print STDERR "$name[0]\n";
        	$name[0] =~s/ID=//;	
        	#print STDERR "$name[0]\n";
        
   	 	if (defined $gene_ID{$name[0]}){
            		for(my $bp = ($gene[3]); $bp < ($gene[4]); $bp++){	 
                		my $position = "$gene[0]_$gene[6]_$bp";
               			$coding{$position} = $name[0];
            		}
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
	open (BAM, "samtools view /path/to/bam/$file\.bam |"); #change file path (files need to be in one folder)
    
    #define outside of the loop and set to 0 (this way the total read number is counted)
    my $read_output = 0;

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
        my $number_of_coding = 0;
        
        #compare the reads from the bam file with the four categories from the gff file (by position and basepair)
        for (my $i = $column[3]; $i < $end; $i++) {
            my $pos_bam = join("_", $column[2], $strand, $i);
            if ($coding{$pos_bam}) {
                $number_of_coding++;
            }
        }

        #only count reads that overlap >= 80% with gff (change percent (e.g. 0.5 for 50%))
        my $math = $number_of_coding / $readlength_clipped;
        if ($math >= 0.8) {$read_output++};
    }
    close BAM;
    # Print the results in output file split by tab (change for "," to get a true .csv file)
    print OUT "$file\t$output\t$read_output\n";
}
close OUT;

exit;
