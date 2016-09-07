#!/usr/bin/perl -w
 
# random_dataset_generate.pl
# This script will create three random datasets as a control to the MAKER GC content project. 

# Megan Bowman
# 28 October 2014

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use POSIX;

my $usage = "\nUsage: $0  --transcript <name of file for est2genome output> --random1 <name of output random file1> --random2 <name of output random file2> --random3 <name of output random file3>\n";

my ($est2genome, $random1, $random2, $random3, $transcript);

GetOptions('transcript=s' => \$transcript, 
	   'random1=s' => \$random1, 
	   'random2=s' => \$random2, 
	   'random3=s' => \$random3);

if (!defined $transcript) { 
    die $usage;
}

if (!-e $transcript) {
  die "$est2genome doesn't exists!\n"
}


my %seq_objects;

my $fasta = Bio::SeqIO->new (-format => 'fasta', -file => $transcript);

my (%gc_hash);

while (my $seqobj2 = $fasta->next_seq()) {
  my $gc_count;
  my $fastaseq = $seqobj2->seq();
  my $seq_id = $seqobj2->display_id();
  my $seq_length = $seqobj2 ->length();
  if ($seq_length >= 800 && $seq_length <= 2400) {
    $seq_objects{$seq_id} = $seqobj2;
    my $gc_percent = rand(100);
    $gc_hash{$seq_id} = $gc_percent;
  }
  next;
}

my $bottom_seqio = Bio::SeqIO->new (-format => 'fasta', -file => ">$random1");
my $top_seqio = Bio::SeqIO->new (-format => 'fasta', -file => ">$random2");
my $med_seqio = Bio::SeqIO->new (-format => 'fasta', -file => ">$random3");



my ($bottom_quartile, $upper_quartile, $middle_quartile, $mid);

my @sorted_ids = sort mysort keys %gc_hash;

my $count = scalar @sorted_ids;
$bottom_quartile = floor ($count/4);
$upper_quartile = $count - $bottom_quartile;
$mid = $bottom_quartile *2;
$middle_quartile = $count - $mid;

for (my $i = 0; $i < $bottom_quartile; ++$i) {
  $bottom_seqio->write_seq($seq_objects{$sorted_ids[$i]});
}

for (my $i = $upper_quartile; $i < $count; ++$i) {
  $top_seqio->write_seq($seq_objects{$sorted_ids[$i]});
}


for (my $i = $bottom_quartile; $i < $middle_quartile; ++$i) {
  $med_seqio->write_seq($seq_objects{$sorted_ids[$i]});
}





exit 0;

sub mysort {
  return $gc_hash{$a} <=> $gc_hash{$b};
}






  
    


      
  
  
  

