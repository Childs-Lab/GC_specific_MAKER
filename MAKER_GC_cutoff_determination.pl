#!/usr/bin/perl -w
 
# MAKER_GC_cutoff_determination.pl
# This script will identify CDS from a genome FASTA and the corresponding MAKER gff, it will then determine the GC content of the fasta and produce two output files, one for making a GC distribution graph and the other cutoffs to be used as input for MAKER_GC_training_set_create.pl 

# Jane Pulman 
# 23 August 2015

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use List::Util qw(sum);
use POSIX;
use Bio::Seq::Quality;
my $usage = "\nUsage: $0  --fasta <full path to file with genome fasta sequences> --gff <full path to maker created GFF3 of est2genome> --name <name to be added to the output files> --peak <peak determination window odd number defaults is 5> --smooth <smoothing window odd number default is  7>\n";

my ($fasta, $gff, $name, $peak, $smooth);

GetOptions('fasta=s' => \$fasta,
	   'gff=s' => \$gff,
	   'name=s' => \$name,
	   'peak=i'=>\$peak,
	   'smooth=i'=>\$smooth);

if (!defined $fasta || !defined $gff) { 
    die $usage;
}

if (!defined $name) {
    $name = "output";
}
if (!defined $peak) {
    $peak = 5;
}
if (!defined $smooth) {
    $smooth = 7;
}
if (!-e $fasta) {
  die "$fasta doesn't exists!\n"
}



if (!-e $gff) {
  die "$gff doesn't exists!\n"
}


#--------------------------------Pull out CDS fasta from gff and fasta input--------------------------------------------------#

my (%fasta_hash, %cds_hash, $new_fasta_header,$check_id , $substring, $strand, $new_id);
my $CDS_OUT = $name."_CDS.fasta";
my $seqin = Bio::SeqIO->new (-file=>"$fasta",-format=>'fasta') ;

my $seqout = Bio::SeqIO->new (-file=>">$CDS_OUT",-format=>'fasta') ;


while (my $seq_obj= $seqin->next_seq) {
    my $seq_id = $seq_obj->id;
    $fasta_hash{$seq_id}=$seq_obj;
}


open my $maker_gff, "$gff" or die " Can't open $gff";
while (my $line = <$maker_gff>) {
    chomp $line;
    if ($line =~ /^##FASTA/ || $line =~ /^>/){
    	# If the gff has the fasta sequence attached to the bottom,
 	# we will crash if we attempt to process the fasta section.
 	last;
    }
    elsif ($line =~ /^#/) {
	next;
    }
    else{
	my @start;
	my @end;
	my @temp_gff = split /\t/, $line;
	if ($temp_gff[2] =~ /mRNA/) {
	    $check_id = $temp_gff[0];
	    $strand = $temp_gff[6];
	    if ($temp_gff[8] =~ /ID=(.+);Parent/) {
		$new_fasta_header = $1;
	    }
	    while (my $line = <$maker_gff>) {
		chomp $line;

		if ($line =~ /^#/){
		    last;
		}

		my @check_cds = split /\t/, $line;
		my $feature = $check_cds[2];

		if ($feature =~/CDS/) {
		    # Now check if this is the first splice form.
 		    # We will not process any splice forms after the first.
 		    # ID=maker-scaffold702-est_gff_est2genome-gene-15.0-mRNA-1:cds;Parent=....
 		    if ($check_cds[8] =~ /mRNA-(\d+):cds;Parent=/) {
 			my $variant_count = $1;

 			print "$check_cds[8]\n";
 			print "$variant_count\n";

 			if ($variant_count != 1) {
 			    # This is not the first splice form.
 			    last;
 			}
 		    }
		    push @start, $check_cds[3];
		    push @end, $check_cds[4];
		}
		elsif (($feature =~ /mRNA/) || ($feature =~ /exon/) || ($feature =~ /contig/) || ($feature =~/three/) || ($feature=~/five/))  {
		    
		}
		else {
		    last;
		}
	    }

	    if (exists $fasta_hash{$check_id}) {

		#######################
		# If this is a positive strand gene, the @start and @end arrays should be sorted from smallest to largest values.
		# If this is a negative strand gene, the @start and @end arrays should be sorted from largest to smallest values.

		my (@sorted_start, @sorted_end);

		if ($strand eq '+') {
		    @sorted_start = sort {$a <=> $b} @start;
		    @sorted_end = sort {$a <=> $b} @end;
		}
		elsif ($strand eq '-') {
		    @sorted_start = sort {$b <=> $a} @start;
		    @sorted_end = sort {$b <=> $a} @end;
		}
		else {
		    die "\nNo strand provided for $new_fasta_header\n\n";
		}

		my $new_fasta;
		for (my $i = 0; $i < @start; $i++) {

		    $substring = $fasta_hash{$check_id}->subseq($sorted_start[$i],$sorted_end[$i]);
		    if ($strand eq '+'){
			$new_fasta .= $substring;
			$new_id = $new_fasta_header."_CDS_+ve";
			next;
		    }
		    elsif ($strand eq '-'){
			my $rev_substring = reverse_comp_subroutine($substring);
			$new_fasta.= $rev_substring;
			$new_id =$new_fasta_header."_CDS_-ve";
			next;
		    }
		}

		my $new_fasta_obj =    Bio::Seq::Quality->new( -id   => $new_id,
							       -seq  => $new_fasta);
		$seqout->write_seq($new_fasta_obj);

	    }
	    else {
		print "$check_id has no sequence in the fasta file\n";
	    }
	}
	
	
    }
}
#--------------------------------Determine GC content from fasta--------------------------------------------------#

my $seqin2 = Bio::SeqIO->new (-file=>"$CDS_OUT",-format=>'fasta') ;
my $GC_OUT = $name."_gc_content.txt";
open my $out, ">$GC_OUT" or die "can't open GC output file" ;
my $total_gc_count = 0;
my $total_A = 0;
my $total_T =  0;
my $total_G =  0;
my $total_C =  0;
my $total_N = 0;
my $total_all = 0;

print $out "seq_id\tAs\tTs\tCs\tGs\ttotal\tGC\%\n";
while (my $seq_obj= $seqin2->next_seq) {
    my $A = 0;
    my $T = 0;
    my $C = 0;
    my $G = 0;
    my $N = 0;
    my $GC = 0;
    my $total = 0;
    my $seq_id = $seq_obj->id;
    my $seq = $seq_obj->seq;
    my @nucs_array = split //, $seq;
    foreach my $nucleotide (@nucs_array){
        if ($nucleotide eq "A"){
            $A++;
        }
        elsif ($nucleotide eq "T"){
            $T++;
	}
        elsif ($nucleotide eq "C"){
            $C++;
        }
        elsif ($nucleotide eq "G"){
            $G++;
        }
        else {
          #  print "found non standard nucleotide sequence $nucleotide \n";
        }

    }


    $total_A += $A;
    $total_T += $T;
    $total_C += $C;
    $total_G += $G;
    $GC += $G + $C;
    $total_gc_count += $GC;
    $total = $A + $C + $G + $T;
    $total_all = $total_all + $A + $C + $G + $T;
    my $GC_content = ($GC/$total)*100;
    my $rounded_GC_content = int($GC_content + 0.5);
    print $out "$seq_id\t$A\t$T\t$C\t$G\t$total\t$GC_content\t$rounded_GC_content\n";

}

#-------------------------------bin GC and determine cut offs for MAKER GC protocol--------------------------------------------------#

my (%hash, %binned_hash, %percent_hash, @values,  $peak_index, $peak2_index, @peaks, @newarray, @sorting);
my $DIST_OUT = $name."_gc_distribution.txt";
my $CUTOFF_OUT = $name."_cutoff.txt";
open my $codonw_in_fh, "$GC_OUT" or die "can't open GC output file for input";
open my $output_fh, ">$DIST_OUT" or die "cant open distribution output";
open my $output2_fh, ">$CUTOFF_OUT" or die "cant open cut off output file";



## this section bins the GC into integers and does a total gene count for each bin

my $total_number_genes = 0;

while (my $gene = <$codonw_in_fh>) {
    ++$total_number_genes;
    chomp $gene;
    if ($gene =~ /GC/) {
        next;
    }
    else {
        my @split_gene = split /\t/, $gene;
        # times the following by 100 for codon and use index 10
	if (!exists($split_gene[6])) {
	    next;
	}
        my $GC = ($split_gene[6]);
	#print "$gene\n";
	#print "$GC\n";
        my $rounded_GC = int ($GC+0.5);
        if (!exists $hash{$rounded_GC}) {
            $hash{$rounded_GC} = 1;
        }
        else {
            $hash{$rounded_GC}++;
        }
    }
}

## this section takes the bins and works out an average looking at bins either side as denoted by user -s input
my @ordered_list = sort keys %hash;

my $last = pop @ordered_list;
my $first = shift @ordered_list;


my $testrange = floor($smooth/2);

for (my $bin = $first; $bin <= $last; ++$bin) {
    my @bin;
    if (exists $hash{$bin}){
        push @bin, $hash{$bin};
    }
    for (my $y= 1; $y <= $testrange; $y++) {

        if (exists $hash{$bin+$y}) {
            push @bin, $hash{$bin+$y};
        }
        if (exists $hash{$bin-$y}) {
            push @bin, $hash{$bin-$y};
        }

    }
    if (@bin) {
        my $average = mean(@bin);
        $binned_hash{$bin}= $average;
    }
}
## this works out a percentage of genes for each bin's average number of genes

print $output_fh "GC\tpercent of total genes\n";
my @output;
foreach my $keys (keys (%binned_hash)) {
    my $percentage = (($binned_hash{$keys}/$total_number_genes)*100);
    $percent_hash{$keys}=$percentage;
}

my $first_value = $percent_hash{$first};
my $last_value = $percent_hash{$last};

## the following line can be used to check that the percentage hash has worked.
#foreach my $percs (keys (%percent_hash)){
# print $percs."\t".$percent_hash{$percs}."\n";
#}

## this sorts the percent hash in order of GC and puts the GCs into an array @keys
my @keys = sort keys %percent_hash;

## this takes the GC values from @ keys and puts there corresponding percentanges in to a @ value array so the index is the same and prints to output for graphical analysis
foreach my $element (@keys){
    push @values, $percent_hash{$element};
}
for (my $i=0; $i <= $#keys; $i ++){
    print $output_fh $keys[$i]."\t".$values[$i]."\n";
}
## Pull out all peaks in the graph looking at for values that are higher than the 3 values either side of the test value
for (my $test_value = 6; $test_value <= @values; $test_value++){

    my $is_peak = 1;
    my $mid_point = floor($peak / 2);
    for (my $x = 0; $x < $peak; ++$x) {
	print "$test_value\t$mid_point\t$x\n";
        if ($x == $mid_point) {
            next;
        }
	elsif (!defined ($values[$test_value - $mid_point]) || !defined ($values[$test_value - $x])){
	    next;
	}
        elsif ($values[$test_value - $mid_point] <= $values[$test_value - $x]) {
            $is_peak = 0;
        }
    }

    if ($is_peak == 1) {
        $peak_index= $test_value - $mid_point;
        push @peaks, $peak_index;
    }

}
## if there are one or two peaks (as expected in for most grass genomes) this section will print out GC cut offs and it will give a warning if there is only one peak
if (@peaks == 0){
    print "no peaks detected\n";
}
elsif (@peaks == 1){
    print "only 1 peak determined the peak value can be see in $CUTOFF_OUT use $DIST_OUT to confirm graphically\n";
    print $output2_fh $keys[$peaks[0]];
}
elsif (@peaks == 2){
    if (($peaks[0] == 0) || ($peaks[1] == 0)){
        print "this is not a bimodal distribution, the peak value can be seen in $CUTOFF_OUT use $DIST_OUT to confirm graphically\n";
        print $output2_fh $keys[$peaks[0]];
    }
    else{
        print "bimodal distribution: see $CUTOFF_OUT for low and high GC cut offs. confirm with graph using $DIST_OUT\n";
        print $output2_fh $keys[$peaks[0]]."\t". $keys[$peaks[1]];
    }
}
elsif (@peaks >=3){

    ## if there are more than 2 peaks (can happen with rough data) this will put the percentages of each peak into an array @sorting, this array will then be sorted numerically and the highest two values are stored as $p1 and $p2 and used to determine their corresponding index in the percentage array (@values) this can then be used to pull out the corresponding GC as the index should be the same in @keys.

    foreach my $sorts (@peaks){
        push @sorting, $values[$sorts];
    }
    @newarray = sort {$a cmp $b} @sorting;
    my $p1 = pop @newarray;
    my $p2 = pop @newarray;
    my $want = $p1;
    my $index = 0;
    ++$index until $values[$index] == $want or $index > $#values;
    my $want2 = $p2;
    my $index2 = 0;
    ++$index2 until $values[$index2] == $want2 or $index2 > $#values;
    print "bimodal distribution: see $CUTOFF_OUT for low and high GC cut offs. confirm with graph using $DIST_OUT\n";
    print $output2_fh  $keys[$index]."\t ".$keys[$index2];
}
#--------------------------------------------------------------------------------------------------------------------#

exit;

sub reverse_comp_subroutine {
    my ($variable) = @_;
    my$revcomp = reverse($variable);
    $revcomp=~tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}
sub mean {
    return sum(@_)/@_;
}
