
BCNAME="test.fa" #don't include .fastq

### in step (2) when identify target seq in reads
#mapped shortest target length cutoff
export targetmappingcutoff=30
#alignment evalu cutoff
export evaluecutoff=0.000000001
#alignment of target missed end length tolerence. 
#if one end of the target missed too long seq, don't consider 
#the alignment. Increase this value can get more reads with target seq
export target_end_error_tolerence=60
#flanking seuqnece length cutoff
export flanking_length_cutoff=20
##################################

### in step (4) when identify insertions with mapped position of flanking seqs

#unqiue alignments or not
export uniquealg=1 #1 is true for unqiue alingment, 0 or other is false

#for multiple flanking alignments, all seq with similar bitscores will be kept, 
#if the new one's bitscore and the largest one ratio in the cutoff range [0, 1].
export bitscoreratiocutoff=0.8


export alignment_length_cutoff=20
export evalue_cutoff=0.00000001
export bin_size=200
#gap between mapped flanking seq to the end of the flanking seq. change this to large number to get more insertion candidates.
export gap_to_target_cutoff=50
#diff between flanking length and aliged length
export diff_alignemnt_total_length=400
#split whole genome into several bins
export genome_bin=100000

#genome blast lib directory and name
#GENOMEBLASTLIB="/work5/jias/nanopore/soybeangenome/Glycine_max.V1.0.dna.toplevel"
GENOMEBLASTLIB="ref"

#gtf file directory and name
#export gtffilename="/work5/jias/nanopore/soybeangenome/Glycine_max.V1.0.37.gtf"
export gtffilename="ref.gtf"
##################################




#echo "(0) convert fastq to fastq"

#perl fastq_to_fasta.pl "$BCNAME".fastq

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

