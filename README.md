# GC-Specific MAKER
###This repo contains the Perl and Python code from "A modified GC-specific MAKER structural genome annotation method reveals improved and novel gene predictions of high and low GC content in Oryza sativa" Bowman et al, 2017.

![GC MAKER](https://github.com/Childs-Lab/GC_specific_MAKER/blob/master/gc_paper_Figure_3.png "GC Specific MAKER")

**Figure 3. Six HMMs MAKER structural annotation method. The center workflow depicts the standard method for training hidden markov models for use in MAKER, while the low GC (top) and high GC (bottom) training methods can be used after creating high and low GC HMM training data sets. After separately training HMMs with the high and low GC training data, all three SNAP HMMs and all three AUGUSTUS HMMs were specified with the maker_opts.ctl file, and MAKER was run to create the six HMMs annotation, which incorporates gene predictions from the standard, high and low GC MAKER runs.**

The MAKER_GC_cutoff_determination.pl script helps to identify the GC values of the peaks in a bimodal grass gene GC content distribution. The script pulls out the CDS fasta sequences for the transcript-based gene predictions from the GFF3 and calculates the GC content for each gene prediction.  The script assigns the gene GC values to integer bins and writes the results to a file. This output file can be used in R to plot the distribution of gene GC content. A fasta file of the CDS sequences and a GC content file (showing nucleotide composition and GC content of each prediction) are also produced. In addition, a text file is created with the high and low peak values of the bimodal gene GC distribution that serve as set points in creating the high and low GC HMM training sets. However, users may pick their own high-GC and low-GC cutoff values, and the gene GC content distribution graph may aid in picking those cutoff values.

```
MAKER_GC_cutoff_determination.pl  --fasta <full path to file with genome fasta sequences> --gff <full path to maker created GFF3 of est2genome> --name <base name for output files> --peak <peak determination window, odd integer, default is 5> --smooth <smoothing window, odd integer, default is 7>
```
The MAKER_GC_training_set_create.py script will create high-GC and low-GC data sets that can be used for training SNAP or AUGUSTUS. The name parameter serves as a prefix for all output files from this script, which include FASTA files of the selected transcripts for HMM training and the corresponding GFF3 file for these transcripts. 

```
MAKER_GC_training_set_create.py --input_file_gff <path to MAKER GFF file> --input_file_GC_content <BASE_NAME_gc_content.txt file> --input_file_GC_cutoff <BASE_NAME_cutoff.txt file> --output_file_low <path to the low GC GFF file> --output_file_high <path to the high GC GFF file> --genome_fasta <path to the genome fasta file>
```
The name parameter serves as a prefix for all output files from this script, which include FASTA files of the selected transcripts for HMM training and the corresponding GFF3 file for these transcripts. 

It should be noted that during the process of SNAP training the total GC content and therefore the number of genes falling below or above the designated cutoffs may be too few for AUGUSTUS training. For GC specific HMM training in O. sativa, SNAP1 transcripts were used for AUGUSTUS training instead of SNAP2 due to this limitation of available transcripts meeting the threshold for AUGUSTUS training. 

###Training AUGUSTUS 

The train_augustus.sh shell script prepares training and testing data sets and makes use of the autoAug.pl training script from AUGUSTUS to create the appropriate HMM files.  

```
train_augustus.sh <path to working directory for training> <path to MAKER gff3 output from initial MAKER run> <species name for AUGUSTUS HMM directory> <path to single fasta file with all transcript assemblies>
```

==========================================
###Creation of randomized MAKER annotations

To assess the impact of GC specific HMM training on the structural annotation of O. sativa, three MAKER annotations were created using HMMs trained with transcripts from the default annotation with randomized GC content. The following perl scripts were used, which create the training GFF3 files based on a random seed instead of percentage GC content. 

The random_dataset_generate.pl script inputs the MAKER standard transcript FASTA and outputs three transcript FASTAs that are used for downstream GFF3 creation and HMM training.

```
random_dataset_generate.pl --transcript < name of file for fasta_merge results from est2genome datastore> --random1 <name of output random file1> --random2 <name of output random file2> --random3 <name of output random file3>
```
The seq_name.pl script is run for each of the three random output FASTA files, which generates a list of MAKER standard transcript names from each transcript FASTA.
```
seq_name.pl --fastafile <path to an output of random_dataset_generate.pl> --output <name of text file with ID names for each FASTA file>
```
Finally, the random_gff3_create.pl script inputs the MAKER standard GFF3 with the genome FASTA included and generates the final randomized GFF3s that are used for SNAP and AUGUSTUS HMM training. 
```
random_gff3_create.pl --align_gff  <path to MAKER GFF3 with FASTA included> --rand_1 <path to random 1 IDs> --rand_2 <path to random 2 IDs>  --rand_3 <path to random 3 IDs>
```
The outputs of these steps are three GFF3 files containing the coordinates of randomly selected gene predictions. Each of these GFF3 files were then used for SNAP or AUGUSTUS training. 

==========================================
###Creation of MAKER standard gene sets

To identify the high-quality gene set, the MAKER accessory scripts gff_merge and fasta_merge, which are included in the MAKER installation, were are used to generate a gff3 file with all gene predictions and evidence data and the transcript and protein fasta files for those predictions.  Pfam domains are identified within the predicted proteins using [hmmscan](https://www.ebi.ac.uk/interpro/search/sequence-search). 

```
hmmscan --domE 1e-5 -E 1e-5 --tblout <MAKER max predictions hmmscan output file> <path to Pfam-A.hmm> <path to predicted protein fasta file>
```
The annotation GFF3 file, the transcript and protein fasta files and the hmmscan results file were used to generate the final high quality MAKER standard gene set. 

```
generate_maker_standard_gene_list.pl  --input_gff <output of gff3_merge> --pfam_results <hmmscan output> --pfam_cutoff  1e-10  --output_file <path to MAKER standard gene list>

get_subset_of_fastas.pl  -l <path to MAKER standard gene list> -f <fasta_merge output transcript/protein fasta>  -o <path MAKER standard transcript/protein fasta>

create_maker_standard_gff.pl  --input_gff <output of gff3_merge> --output_gff <path to MAKER standard GFF3> --maker_standard_gene_list <path to MAKER standard gene list>
```
==========================================
###Filtering gene sets for transposable elements

Predicted proteins were compared to a database of [Gypsy transposable elements (3.1.b2)](http://gydb.org/index.php/Collection_alignments).  Predicted proteins were also aligned with blastp to a database of [transposases]( http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced). A GFF3 file of TE-related genes was derived from the MSU-RGAP gene annotation GFF3 [file] (http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3) and was compared to the MAKER standard GFF3 file using [gffcompare](https://github.com/gpertea/gffcompare).

```
hmmscan --tblout <Gypsy HMM analysis output> -E 1e-5 --domE 1e-5 <path to gypsy_db_3.1b2.hmm> <path to maker standard proteins fasta> 

blastp -db <Tpases020812 database> -query <path to MAKER standard protein fasta> -out <path to Tpases blast output> -evalue 1e-10 -outfmt 6 

gffcompare -o <TE comparison output file> -r <MSU RGAP TE GFF3> <MAKER standard GFF3>
```
The create_no_TE_genelist.py script use the data derived above, the Pfam hmmscan results file and a list of TE-related Pfam domains (TE_Pfam_domains.txt) to create a list of MAKER standard genes with no TE-related predictions.

```
create_no_TE_genelist.py --input_file_TEpfam <TE_Pfam_domains.txt> --input_file_maxPfam <MAKER max predictions hmmscan output file> --input_file_geneList_toKeep <path to MAKER standard gene list> --input_file_TEhmm <Gypsy HMM analysis output> --input_file_TEblast <path to Tpases blast output> --input_file_TErefmap <TE comparison output refmap file> --output_file <path to TE filtered MAKER standard gene list>
```














