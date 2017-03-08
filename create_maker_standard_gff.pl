#! /usr/bin/perl

# create_maker_standard_gff.pl

# 16 October 2013
# Kevin Childs

###

use Getopt::Long;
use strict;
use warnings;

my $usage = "\n$0\n    --input_gff   <path to input gff file>\n" .
            "    --output_gff  <path to output gff file>\n" .
            "    --maker_standard_gene_list   <path to list of maker standard genes>\n" .
            "    [--help]\n\n";

my ($input_gff, $output_gff, $maker_standard_gene_file);
my $help;

Getopt::Long::GetOptions( "input_gff=s" => \$input_gff,
                          "output_gff=s" => \$output_gff,
                          "maker_standard_gene_list=s" => \$maker_standard_gene_file,
                          "help" => \$help) || die;

if (defined($help)) {
    print $usage;
    exit;
}
if (!defined($input_gff) ||  !(-e $input_gff) ||
    !defined($output_gff) ||  (-e $output_gff) ||
    !defined($maker_standard_gene_file) ||  !(-e $maker_standard_gene_file)) {
    die $usage;
}

my ($gff_in_fh, $gff_out_fh, $gene_list_fh);

open $gff_in_fh, $input_gff or die "Unable to open $input_gff for reading.\n";
open $gff_out_fh, ">$output_gff" or die "Unable to open $output_gff for writing.\n";
open $gene_list_fh, $maker_standard_gene_file or die "Unable to open $maker_standard_gene_file for reading.\n";

my %maker_standard_genes;

while (my $line = <$gene_list_fh>) {
    chomp $line;

    $maker_standard_genes{$line} = 1;

    $line =~ s/-mRNA-\d+$//;
    $line =~ s/-mRNA-\d+-LOW$/-LOW/;
    $line =~ s/-mRNA-\d+-HIGH$/-HIGH/;
    $maker_standard_genes{$line} = 1;

}
close $gene_list_fh;

#3070882	maker	gene	2053	3823	.	+	.	ID=snap_masked-3070882-abinit-gene-0.3;Name=snap_masked-3070882-abinit-gene-0.3;
#3070882	maker	mRNA	2053	3823	.	+	.	ID=snap_masked-3070882-abinit-gene-0.3-mRNA-1;Parent=snap_masked-3070882-abinit-gene-0.3;Name=snap_masked-3070882-abinit-gene-0.3-mRNA-1;_AED=1.00;_eAED=1.00;_QI=0|0|0|0|1|1|3|0|189;
#3070882	maker	exon	2053	2120	-9.442	+	.	ID=snap_masked-3070882-abinit-gene-0.3-mRNA-1:exon:564177;Parent=snap_masked-3070882-abinit-gene-0.3-mRNA-1;
#3070882	maker	exon	2422	2595	-9.104	+	.	ID=snap_masked-3070882-abinit-gene-0.3-mRNA-1:exon:564178;Parent=snap_masked-3070882-abinit-gene-0.3-mRNA-1;
#3070882	maker	exon	3496	3823	17.396	+	.	ID=snap_masked-3070882-abinit-gene-0.3-mRNA-1:exon:564179;Parent=snap_masked-3070882-abinit-gene-0.3-mRNA-1;
#3070882	maker	CDS	2053	2120	.	+	0	ID=snap_masked-3070882-abinit-gene-0.3-mRNA-1:cds:557580;Parent=snap_masked-3070882-abinit-gene-0.3-mRNA-1;
#3070882	maker	CDS	2422	2595	.	+	1	ID=snap_masked-3070882-abinit-gene-0.3-mRNA-1:cds:557581;Parent=snap_masked-3070882-abinit-gene-0.3-mRNA-1;
#3070882	maker	CDS	3496	3823	.	+	1	ID=snap_masked-3070882-abinit-gene-0.3-mRNA-1:cds:557582;Parent=snap_masked-3070882-abinit-gene-0.3-mRNA-1;

while (my $line = <$gff_in_fh>) {
    chomp $line;

    if ($line =~ /^##FASTA/) {
	last;
    }
    if ($line =~ /^#/) {
	print $gff_out_fh "$line\n";
	next;
    }

    my @elems = split "\t", $line;

    if ($elems[1] ne 'maker') {
	# All non-maker features will be retained.
	print $gff_out_fh "$line\n";
	next;
    }

    if ($elems[2] eq 'gene') {
	if ($elems[8] =~ /Name=(.+)/) {
	    my $id = $1;

	    if (exists($maker_standard_genes{$id})) {
		print $gff_out_fh "$line\n";
	    }
	}
	next;
    }

    if ($elems[2] eq 'mRNA') {
	if ($elems[8] =~ /Parent=(.+);Name=(.+);/) {
	    my $parent_id = $1;
	    my $id = $2;

	    if (exists($maker_standard_genes{$parent_id})) {
		print $gff_out_fh "$line\n";
		$maker_standard_genes{$id} = 1;
	    }
	}
	next;
    }

    if ($elems[2] eq 'exon' || $elems[2] eq 'CDS' || $elems[2] =~ /UTR/) {
	if ($elems[8] =~ /Parent=(.+)/) {
	    my $parent_id = $1;

	    if (exists($maker_standard_genes{$parent_id})) {
		print $gff_out_fh "$line\n";
	    }
	}
	next;
    }

}
close $gff_in_fh;
close $gff_out_fh;

exit;
