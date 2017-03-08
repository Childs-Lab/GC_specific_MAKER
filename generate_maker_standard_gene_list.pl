#! /usr/bin/perl

# generate_maker_standard_gene_list.pl

# 20 September 2013
# Kevin Childs

use Getopt::Long;
use strict;
use warnings;

# We ain't waiting on no stinking print buffer.
$| = 1;

my $usage = "\n$0\n    --input_gff   <path_to_input_gff_file>\n" .
            "    --pfam_results  <path_to_pfam_results_file>\n" .
            "    --pfam_cutoff   <pfam p-value cutoff>\n" .
            "    --output_file   <MAKER-standard gene list output file>\n" .
            "    [--help]\n\n";

my ($input_gff, $pfam_results, $pfam_cutoff, $output_file);
my $help;

Getopt::Long::GetOptions( "input_gff=s" => \$input_gff,
                          "pfam_results=s" => \$pfam_results,
                          "pfam_cutoff=f" => \$pfam_cutoff,
                          "output_file=s" => \$output_file,
                          "help" => \$help) || die;

if (defined($help)) {
    print $usage;
    exit;
}
if (!defined($input_gff) ||  !(-e $input_gff) ||
    !defined($pfam_results) ||  !(-e $pfam_results) ||
    !defined($pfam_cutoff) || 
    !defined($output_file) ||  (-e $output_file)) {
    die $usage;
}

my @bad_terms;
#my @bad_terms = ("Retrotrans_gag",
		 #"RVT_",
		 #"Transposase",
		 #"gag_pre-integrs",
		 #"rve",
		 #"UBN2",
		 #"RNase_H",
    #);

open PFAM, $pfam_results or die "\nUnable to open $pfam_results for reading.\n\n";

my %pfam_genes;
my %bad_genes;
while (my $line = <PFAM>) {
    chomp $line;

    if ($line =~ /^#/) {
	next;
    }

    # There are no tabs in these result lines.
    # Put the tabs in there.
    $line =~ s/\s+/\t/g;

    my @elems = split "\t", $line;

    my $is_bad_domain = 0;
    foreach my $domain_name (@bad_terms) {
	if ($elems[0] =~ /^$domain_name/) {
	    $is_bad_domain = 1;
	    last;
	}
    }
    if ($is_bad_domain == 1) {
	$bad_genes{$elems[2]} = 1;
	next;
    }

    #print "e-value $elems[4]\n";
    if ($elems[4] <= $pfam_cutoff) {
	$pfam_genes{$elems[2]} = $elems[4];
	#print "$elems[2]\n";
    }
}
close PFAM;

#print "number of pfam genes: " . scalar(keys(%pfam_genes)) . "\n";

open GFF, $input_gff or die "\nUnable to open $input_gff for reading.\n\n";
open OUT, ">$output_file" or die "\nUnable to open $output_file for writing.\n\n";

my %maker_standard_genes;
while (my $line = <GFF>) {
    chomp $line;

    if ($line =~ /#FASTA/) {
	last;
    }
    if ($line =~ /^#/) {
	next;
    }

    my @elems = split "\t", $line;

    if ($elems[1] eq 'maker' && $elems[2] eq 'mRNA') {
	if ($elems[8] =~ /^ID=([^;]+)/) {
            my $id = $1;
            if ($elems[8] =~ /_AED=([\.\d]+);/) {
		my $aed = $1;
		#print "$aed\t$id\n";
		if (!exists($bad_genes{$id}) && ($aed < 1 || exists($pfam_genes{$id}))) {
		    # This is the MAKER-standard gene set.
		    print OUT "$id\n";
		}
            }
            else {
		print "Failed to regex AED value.\n$line\n";
            }
	}
    }

}
close GFF;
close OUT;

exit;
