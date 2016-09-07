#!/usr/bin/perl -w
 
# random_gff3_create.pl
# This script will create new gff3s for training SNAP based on random datasets.

# Megan Bowman
# 29 October 2014

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $usage = "\nUsage: $0  --align_gff <path to the MAKER alignment gff3> --rand_1 <path to the random 1 IDs file> --rand_2 <path to the random 2 IDs file> --middle_quart <path to the random 3 IDs file>\n";

my ($align_gff, $rand_1, $rand_2, $rand_3);

GetOptions('align_gff=s' => \$align_gff,
	   'rand_1=s' => \$rand_1,
	   'rand_2=s' => \$rand_2,
	   'rand_3=s' => \$rand_3);

if ((!defined $align_gff) || (!defined $rand_1) || (!defined $rand_2) || (!defined $rand_3)) { 
    die $usage;
}

if ((!-e $align_gff) || (!-e $rand_1) || (!-e $rand_2) || (!-e $rand_3)) {
  die $usage;
}

open IN1, "$rand_1" or die "Cannot open rand_1 ids for reading\n";
open IN2, "$rand_2" or die "Cannot open rand_2 ids for reading\n";
open IN4, "$rand_3" or die "Cannot open rand_3 ids for reading\n";


my ($line, %random_1_gene_id_hash, %random_2_gene_id_hash, %random_3_gene_id_hash, %random_1_id_hash, %random_2_id_hash, %random_3_id_hash);

while ($line = <IN1>) {
  chomp $line;
  if ($line =~ />(.+)/) {
    my $id = $1;
    $random_1_id_hash{$id} = 1;  
  }
  if ($line =~ /^>(.+)\-mRNA\-\d$/) {
    my $gene_id = $1;
    $random_1_gene_id_hash{$gene_id} =1;
  }
}

while ($line = <IN2>) {
  chomp $line;
  if ($line =~ />(.+)/) {
    my $id = $1;
    $random_2_id_hash{$id} = 1;
  }
  if ($line =~ /^>(.+)\-mRNA\-\d$/) {
    my $gene_id = $1;
    $random_2_gene_id_hash{$gene_id} =1;
  }
}


while ($line = <IN4>) {
  chomp $line;
  if ($line =~ />(.+)/) {
    my $id = $1;
    $random_3_id_hash{$id} = 1;
  }
  if ($line =~ /^>(.+)\-mRNA\-\d$/) {
    my $gene_id = $1;
    $random_3_gene_id_hash{$gene_id} =1;
  }
}



open IN3, "$align_gff" or die "Cannot open alignment gff3 for reading\n";
open OUT1, ">random_1.gff3" or die "Cannot open random_1 gff3 for writing\n";
open OUT2, ">random_2.gff3" or die "Cannot open random_2 gff3 for writing\n";
open OUT3, ">random_3.gff3" or die "Cannot open random_3 gff3 for writing\n";

print OUT1 "##gff-version 3\n";
print OUT2 "##gff-version 3\n";
print OUT3 "##gff-version 3\n";

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
	if (exists $random_1_gene_id_hash{$gene_id}) {
	  print OUT1 $contig;
	  print OUT1 $line;
	  $contig = "";
	  next;
	}
	if (exists $random_1_id_hash{$gene_id}) {
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
	if (exists $random_1_gene_id_hash{$gene_id}) {
	  print OUT1 $line;
	  next;
	}
	if (exists $random_1_id_hash{$gene_id}) {
	  print OUT1 $line;
	  next;
	}
      }
    }
    if (($elems[2] =~ /exon/) || ($elems[2] =~ /CDS/) || ($elems[2] =~ /five_prime_UTR/) || ($elems[2] =~ /three_prime_UTR/)) {
      if (($elems[8] =~ /Parent=(.+)$/) || ($elems[8] =~ /Parent=(.+);/) || ($elems[8] =~ /ID=(.+);/)) { 
	my $parent_id = $1;
	if (exists $random_1_id_hash{$parent_id}) {
	  print OUT1 $line;
	}
	if ($parent_id =~ /,/) {
	  my @id_array = split ",", $parent_id;
	  for (my $i = 0; $i < $#id_array; $i++) {
	    if (exists $random_1_id_hash{$id_array[$i]}) {
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
	if (exists $random_2_gene_id_hash{$gene_id}) {
	  print OUT2 $contig2;
	  print OUT2 $line2;
	  $contig2 = "";
	  next;
	}
	if (exists $random_2_id_hash{$gene_id}) {
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
	if (exists $random_2_gene_id_hash{$gene_id}) {
	  print OUT2 $line2;
	  next;
	}
	if (exists $random_2_id_hash{$gene_id}) {
	  print OUT2 $line2;
	  next;
	}
      }
    }
    if (($elems2[2] =~ /exon/) || ($elems2[2] =~ /CDS/) || ($elems2[2] =~ /five_prime_UTR/) || ($elems2[2] =~ /three_prime_UTR/)) {
      if (($elems2[8] =~ /Parent=(.+)$/) || ($elems2[8] =~ /Parent=(.+);/) || ($elems2[8] =~ /ID=(.+);/)) {
	my $parent_id = $1;
	if (exists $random_2_id_hash{$parent_id}) {
	  print OUT2 $line2;
	}
	if ($parent_id =~ /,/) {
	  my @id_array2 = split ",", $parent_id;
	  for (my $i = 0; $i < $#id_array2; $i++) {
	    if (exists $random_2_id_hash{$id_array2[$i]}) {
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

close IN3;

open IN3, "$align_gff" or die "Cannot open alignment gff3 for reading\n";

my ($contig3, @elems3);

while (my $line2 = <IN3>) {
  if ($line2 =~ /^\#\#FASTA$/) {
    print OUT3 $line2;
    while (my $line2 = <IN3>) {
      print OUT3 $line2;
    }
  }
  if ($line2 =~ /^#/) {
    next;
  }
  @elems3 = split "\t", $line2;
  if ($elems3[2] =~ /contig/) {
    $contig3 = $line2;
    next;
  }
  if ($elems3[1] =~ /maker/) {
    if ($elems3[2] =~ /gene/) {
      if ($elems3[8] =~ /^ID=(.+);/) {
	my $gene_id = $1;
	if (exists $random_3_gene_id_hash{$gene_id}) {
	  print OUT3 $contig2;
	  print OUT3 $line2;
	  $contig3 = "";
	  next;
	}
	if (exists $random_3_id_hash{$gene_id}) {
	  print OUT3 $contig2;
	  print OUT3 $line2;
	  $contig3 = "";
	  next;
	}
      }
    }
    if ($elems3[2] =~ /mRNA/) {
      if ($elems3[8] =~ /^ID=(.+);Parent=/) {
	my $gene_id = $1;
	if (exists $random_3_gene_id_hash{$gene_id}) {
	  print OUT3 $line2;
	  next;
	}
	if (exists $random_3_id_hash{$gene_id}) {
	  print OUT3 $line2;
	  next;
	}
      }
    }
    if (($elems3[2] =~ /exon/) || ($elems3[2] =~ /CDS/) || ($elems3[2] =~ /five_prime_UTR/) || ($elems3[2] =~ /three_prime_UTR/)) {
      if (($elems3[8] =~ /Parent=(.+)$/) || ($elems3[8] =~ /Parent=(.+);/) || ($elems3[8] =~ /ID=(.+);/)) { 
	my $parent_id = $1;
	if (exists $random_3_id_hash{$parent_id}) {
	  print OUT3 $line2;
	}
	if ($parent_id =~ /,/) {
	  my @id_array2 = split ",", $parent_id;
	  for (my $i = 0; $i < $#id_array2; $i++) {
	    if (exists $random_3_id_hash{$id_array2[$i]}) {
	      print OUT3 $line2;
	      next;
	    }
	    next;
	  }
	}
      }
    }
  }
}  
