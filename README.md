# GC Specific MAKER
*This repo contains the Perl code from "A modified GC-specific MAKER structural genome annotation method reveals improved and novel gene predictions of high and low GC content in Oryza sativa" Bowman et al, 2016.*

The MAKER_GC_cutoff_determination.pl script helps to identify the GC values of the peaks in a bimodal grass gene GC content distribution. The script pulls out the CDS fasta sequences for the transcript-based gene predictions from the GFF3 and calculates the GC content for each gene prediction.  The script assigns the gene GC values to integer bins and writes the results to a file. This output file can be used in R to plot the distribution of gene GC content. A fasta file of the CDS sequences and a GC content file (showing nucleotide composition and GC content of each prediction) are also produced. In addition, a text file is created with the high and low peak values of the bimodal gene GC distribution that serve as set points in creating the high and low GC HMM training sets. However, users may pick their own high-GC and low-GC cutoff values, and the gene GC content distribution graph may aid in picking those cutoff values.
```
MAKER_GC_cutoff_determination.pl  --fasta <full path to file with genome fasta sequences> --gff <full path to maker created GFF3 of est2genome> --name <base name for output files> --peak <peak determination window, odd integer, default is 5> --smooth <smoothing window, odd integer, default is 7>
```
The MAKER_GC_training_set_create.pl script will create high-GC and low-GC data sets that can be used for training SNAP or AUGUSTUS. The name parameter serves as a prefix for all output files from this script, which include FASTA files of the selected transcripts for HMM training and the corresponding GFF3 file for these transcripts. 

```
MAKER_GC_training_set_create.pl  --transcript <path to predicted transcripts fasta file>  --align_gff <full path to maker created GFF3 of est2genome with FASTA included> --name <prefix for the output files> --cutoffs <path to the cutoff output file>
```

The name parameter serves as a prefix for all output files from this script, which include FASTA files of the selected transcripts for HMM training and the corresponding GFF3 file for these transcripts. 

It should be noted that during the process of SNAP training the total GC content and therefore the number of genes falling below or above the designated cutoffs may be too few for AUGUSTUS training. For GC specific HMM training in O. sativa, SNAP1 transcripts were used for AUGUSTUS training instead of SNAP2 due to this limitation of available transcripts meeting the threshold for AUGUSTUS training. 
