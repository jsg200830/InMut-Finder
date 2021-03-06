# Insertion-in-mutagenesis finder, InMut-Finder

# Introduction
Transposon induced mutagenesis technology has been widely used in functional genomics, especially in reverse genetics. The low sequencing cost make transposon mapping realistic. However, the precision of transposon identification is limited due to read length and coverage depth in the next-generation sequencing technology. To overcome this problem, we develop a convenient and fast bioinformatic software tool, named as insertion-in-mutagenesis finder (InMut-Finder), to map the insertion site of the transposons based on the flanking sequences. It is designed to identify both the transposon and flanking sequences from the third sequencing reads (such as Nanopore), and determine the relative genomic coordinates of the transposons and the neighbor genes. It is also allowed to be used to identify the insertion sites of other transgenic elements, for example, T-DNA and residue sequences in transgenic materials, only when the inserted sequences were known. 

In addition, InMut-Finder is designed for the long reads by Nanopore or PacBio, and might be used for HiSeq short reads, as it has been tested successfull in a line of Medicago truncatula with Tnt1 insertion. 

# Dependencies required
InMut-Finder could be run in any of the systems, including windows, Mac and Linux, after the following dependencies are installed. 
        1.Perl
        2.R
        3.BlastN

# Running procedures
The shell script of “run_command.sh” integrated all the parameters and running commands. Simply type “./run_command.sh” in the command line to run all the five steps, after all the input files and parameters are set up in this shell script.

        ./run_command.sh

It will call “identify_target_in_reads_uniqreadid.pl” to search for the long reads covering both insertion target fragment and flanking sequences. The file of “identify_flanking_in_genome_uniq.pl” works on screening the whole genome for the genomic coordinates of insertion and neighboring genes. The R code of “cal_pscore.r” calculates the P values for each insertion, and outputs the final file of “test.fa_insertion_list_nu.txt_score.csv”. 

There are two Perl scripts, one R script and one shell script (run_command.sh), to do the five steps. 

        echo "(1) blast target in reads"
        blastn -query target.fa -subject "$BCNAME"  -out "$BCNAME".blast_table  -max_target_seqs 1000000 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore'
        echo "(2) find targe and get flankning seqs"
        perl identify_target_in_reads_uniqreadid.pl "$BCNAME".blast_table "$BCNAME" > "$BCNAME"_flanking.fasta
        echo "(3) blast flankning seqs"
        blastn -db "$GENOMEBLASTLIB"  -query "$BCNAME"_flanking.fasta  -out "$BCNAME"_flanking.blast_table  -num_threads 30 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send evalue bitscore'   -max_target_seqs 1000000
        echo "(4) find insertions "
        perl identify_flanking_in_genome_uniq.pl "$BCNAME"_flanking.blast_table "$BCNAME"_flanking.fasta "$BCNAME"_insertion_flanking_seqs_nu.txt  "$BCNAME"_insertion_list_nu.txt
        echo "(5) calculate scores"
        Rscript --vanilla cal_pscore.r "$BCNAME"_insertion_list_nu.txt



# Input
The four input files are required to run InMut-Finder in the command line. 

The first one is the fasta file for the inserted target element, for example, Tnt1 or Ds transposons (see the example file of “target.fa”). 

The second is the fasta file which contains all the long reads produced by Nanopore and PacBio sequencing. They might be the re-sequencing of the whole genomes, or the enriched genome sequencing data just like that in the published paper (Li et al., 2020). See the example file of “test.fa”. 

The third is the reference genome, with a fasta format, and the fourth is the genomic annotation file, gtf or gff format. See the example files of “ref.fa” and “ref.gtf”.

# Output
InMut-Finder generated three output files. 

The first one contains the summary of insertion sites, with five columns, which indicate the information of location (genomic coordinates of insertion), gap (extend of insertion), counts (read counts to support this insertion), gene (neighbor genes within a distance of 2000 bp) and P-value (probability value to reject the hypothesis in Poisson distribution). See the example file of “test.fa_insertion_list_nu.txt_score.csv”.

The second contains the 5’ and 3’ sequences from the Nanopore/PacBio sequencing for each insertion. The sequences could be used for designing primers to validate the insertion. See the example file of “test.fa_insertion_flanking_seqs_nu.txt”.

The third contains all the reads which cover both the inserted target fragment and the flanking sequences which can be mapped to the reference genome. This file is useful to do the alignment against the reference genome, for example by minimap2. The bam file can be viewed in IGV tool for manually examining the insertion sites. See the example file of “test.fa_flanking.fasta”. 

# About
InMut-Finder is developped by both Shangang Jia in China Agricultural University and Chi Zhang in University of Nebraska - Lincoln. 

# Citing InMut-Finder
Rui Song, Ziyao Wang, Hui Wang, Han Zhang, Xuemeng Wang, Hanh Nguyen, David Holding, Bin Yu, Tom Clemente, Shangang Jia, Chi Zhang. InMut-Finder: a software tool for insertion identification in mutagenesis using Nanopore long reads. BMC Genomics, 2021, 22: 908. 
