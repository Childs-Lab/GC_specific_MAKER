#!/usr/bin/perl -w
 
# MAKER_GC_training_set_create.pl
# This script will identify est2genome aligned regions in a genome and calculate the GC content. Two fastas are created based on the peak GC content as determined by the GC distribution. FASTAs and GFF3s are created for the high and low GC training sets.  

# Megan Bowman
# 20 March 2015, 19 August 2015

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use POSIX;

my $usage = "\nUsage: $0  --transcript <name of file for est2genome output> --align_gff <full path to maker created GFF3 of est2genome> --name <name to be added to the output files> --cutoffs <path to the cutoff output file>\n";

my ($transcript, $cutoff, $align_gff, $name);

GetOptions('transcript=s' => \$transcript,
	   'cutoffs=s' => \$cutoff,
	   'align_gff=s' => \$align_gff,
	   'name=s' => \$name);

if (!defined $transcript || !defined $cutoff || !defined $align_gff) { 
    die $usage;
}

if (!defined $name) {
    $name = "output";
}

if (!-e $transcript) {
  die "$transcript doesn't exists!\n"
}

if (!-e $cutoff) {
  die "$cutoff doesn't exists!\n"
}

if (!-e $align_gff) {
  die "$align_gff doesn't exists!\n"
}


#--------------------------------Calculate GC Content of est2genome transcripts--------------------------------------------------#

my (%seq_objects, %gc_hash, @cutoffs, @gcs, $filename1, $filename2);

my $fasta = Bio::SeqIO->new (-format => 'fasta', -file => $transcript);

while (my $seqobj2 = $fasta->next_seq()) {
  my $gc_count;
  my $fastaseq = $seqobj2->seq();
  my $seq_id = $seqobj2->display_id();
  my $seq_length = $seqobj2 ->length();
  if ($seq_length >= 800 && $seq_length <= 2400) {
    $seq_objects{$seq_id} = $seqobj2;
    while ($fastaseq =~ /([G|C]+)/g) {
      $gc_count += length($1);
    }
    my $gc_percent = floor(($gc_count/$seq_length) * 100);
    $gc_hash{$seq_id} = $gc_percent;
    push (@gcs, $gc_percent);
  }
  next;
}

#------------------Detect number of peaks in GC distribution and create FASTAs of GC training sets-------------------------------#

$filename1 = $name . "_low_GC_training_set.fasta";
$filename2 = $name . "_high_GC_training_set.fasta";

my $bottom_seqio = Bio::SeqIO->new (-format => 'fasta', -file => ">$filename1");
my $top_seqio = Bio::SeqIO->new (-format => 'fasta', -file => ">$filename2");

open IN, "$cutoff" or die "\nFailed to load file of cutoffs: $!\n";

my ($high, $low, $peak);

while (my $line = <IN>) {
  chomp $line;
  @cutoffs = split "\t", $line;
  if (scalar @cutoffs > 1) {
    print "Detecting bimodal distribution peaks\n";
    $low = $cutoffs[0];
    $high = $cutoffs[1];
    foreach my $key (keys %gc_hash) {
      if ($gc_hash{$key} <= $low) {
	$bottom_seqio->write_seq($seq_objects{$key});
      }
      if ($gc_hash{$key} >= $high) {
	$top_seqio->write_seq($seq_objects{$key});
      }
    }
  }
  if (scalar @cutoffs == 1) {
    print "Detecting a single distribution peak\n";
    $peak = $cutoffs[0];
    foreach my $key (keys %gc_hash) {
      if ($gc_hash{$key} <= $low) {
	$bottom_seqio->write_seq($seq_objects{$key});
      }
      if ($gc_hash{$key} >= $high) {
	$top_seqio->write_seq($seq_objects{$key});
      }
    }
  }
}


#------------------------------------------Creating GFF3s of GC HMM training --------------------------------------------------#

print "Creating high and low GC training GFF3s\n";

my (%top_id_hash_gene, %top_id_hash_mRNA, %low_id_hash_gene, %low_id_hash_mRNA, $low_gff, $high_gff);


my $low_fasta = Bio::SeqIO->new (-format => 'fasta', -file => $filename1);
my $high_fasta = Bio::SeqIO->new (-format => 'fasta', -file => $filename2);

while (my $seqobj = $low_fasta ->next_seq()) {
    my $low_id = $seqobj->display_id();
    $low_id_hash_mRNA{$low_id} = 1;
    $low_id =~ /(.+)\-mRNA\-1/;
    my $gene_id = $1;
    $low_id_hash_gene{$gene_id} = 1;
}


while (my $seqobj = $high_fasta ->next_seq()) {
    my $top_id = $seqobj->display_id();
    $top_id_hash_mRNA{$top_id} = 1;
    $top_id =~ /(.+)\-mRNA\-1/;
    my $gene_id = $1;
    $top_id_hash_gene{$gene_id} = 1;
}

$low_gff = $name . "_low_GC_training_set.gff3";
$high_gff = $name . "_high_GC_training_set.gff3";

open IN3, "$align_gff" or die "Cannot open alignment gff3 for reading\n";
open OUT1, ">$low_gff" or die "Cannot open low GC gff3 for writing\n";
open OUT2, ">$high_gff" or die "Cannot open high GC gff3 for writing\n";

print OUT1 "##gff-version 3\n";
print OUT2 "##gff-version 3\n";

my ($contig, @elems);

while (my $line = <IN3>) {
  if ($line =~ /^\#\#FASTA$/) {
    print OUT1 $line;
    while (my $line = <IN3>) {
      print OUT1 $line;
    }
  }
  if ($line =~ /^#/) {
    next;
  }
  @elems = split "\t", $line;
  if ($elems[2] =~ /contig/) {
    $contig = $line;
    next;
  }
  if ($elems[1] =~ /maker/) {
    if ($elems[2] =~ /gene/) {
      if ($elems[8] =~ /^ID=(.+);/) {
        my $gene_id = $1;
        if (exists $low_id_hash_gene{$gene_id}) {
          print OUT1 $contig;
          print OUT1 $line;
          $contig = "";
          next;
        }
      }
    }
    if ($elems[2] =~ /mRNA/) {
      if ($elems[8] =~ /^ID=(.+);Parent=/) {
        my $gene_id = $1;
        if (exists $low_id_hash_mRNA{$gene_id}) {
          print OUT1 $line;
          next;
        }
      }
    }
    if (($elems[2] =~ /exon/) || ($elems[2] =~ /CDS/) || ($elems[2] =~ /five_prime_UTR/) || ($elems[2] =~ /three_prime_UTR/)) {
      if (($elems[8] =~ /Parent=(.+)$/) || ($elems[8] =~ /Parent=(.+);/) || ($elems[8] =~ /ID=(.+);/)) {
        my $parent_id = $1;
        if (exists $low_id_hash_mRNA{$parent_id}) {
          print OUT1 $line;
        }
        if ($parent_id =~ /,/) {
          my @id_array = split ",", $parent_id;
          for (my $i = 0; $i < $#id_array; $i++) {
            if (exists $low_id_hash_mRNA{$id_array[$i]}) {
              print OUT1 $line;
              next;
            }
            next;
          }
        }
      }
    }
  }
}


close IN3;

open IN3, "$align_gff" or die "Cannot open alignment gff3 for reading\n";

my ($contig2, @elems2);

while (my $line2 = <IN3>) {
  if ($line2 =~ /^\#\#FASTA$/) {
    print OUT2 $line2;
    while (my $line2 = <IN3>) {
      print OUT2 $line2;
    }
  }
  if ($line2 =~ /^#/) {
    next;
  }
  @elems2 = split "\t", $line2;
  if ($elems2[2] =~ /contig/) {
    $contig2 = $line2;
    next;
  }
  if ($elems2[1] =~ /maker/) {
    if ($elems2[2] =~ /gene/) {
      if ($elems2[8] =~ /^ID=(.+);/) {
        my $gene_id = $1;
        if (exists $top_id_hash_gene{$gene_id}) {
          print OUT2 $contig2;
          print OUT2 $line2;
          $contig2 = "";
          next;
        }
      }
    }
    if ($elems2[2] =~ /mRNA/) {
      if ($elems2[8] =~ /^ID=(.+);Parent=/) {
        my $gene_id = $1;
        if (exists $top_id_hash_mRNA{$gene_id}) {
          print OUT2 $line2;
          next;
        }
      }
    }
    if (($elems2[2] =~ /exon/) || ($elems2[2] =~ /CDS/) || ($elems2[2] =~ /five_prime_UTR/) || ($elems2[2] =~ /three_prime_UTR/)) {
      if (($elems2[8] =~ /Parent=(.+)$/) || ($elems2[8] =~ /Parent=(.+);/) || ($elems2[8] =~ /ID=(.+);/)) {
        my $parent_id = $1;
        if (exists $top_id_hash_mRNA{$parent_id}) {
          print OUT2 $line2;
        }
        if ($parent_id =~ /,/) {
          my @id_array2 = split ",", $parent_id;
          for (my $i = 0; $i < $#id_array2; $i++) {
            if (exists $top_id_hash_mRNA{$id_array2[$i]}) {
              print OUT2 $line2;
              next;
            }
            next;
          }
        }
      }
    }
  }
}


  
  

