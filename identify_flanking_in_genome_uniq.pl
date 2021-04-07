#!/usr/bin/perl

# input one blast file against target.fa, and further blast against soybean genome for coordinates and inserted genes.

use warnings;

#unqiue alignments or not
$uniquealg=$ENV{'uniquealg'}; #1 is true, 0 or other is false

#for multiple flanking alignments, all seq with similar bitscores will be kept, 
#if the new one's bitscore and the largest one ratio in the cutoff range
$bitscoreratiocutoff = $ENV{'bitscoreratiocutoff'}; 

#$alignment_length_cutoff = 20;
$alignment_length_cutoff = $ENV{'alignment_length_cutoff'};
#$evalue_cutoff=0.00000001;
$evalue_cutoff= $ENV{'evalue_cutoff'};
#$bin_size=200;
$bin_size=$ENV{'bin_size'};   #insertion location cutting edge boundary bin
#$gap_to_target_cutoff=30;
$gap_to_target_cutoff=$ENV{'gap_to_target_cutoff'}; #gap between mapped flanking seq to the end of the flanking seq.
#$diff_alignemnt_total_length = 200;
$diff_alignemnt_total_length = $ENV{'diff_alignemnt_total_length'}; #diff between flanking length and aliged length
#$genome_bin=100000;
$genome_bin=$ENV{'genome_bin'};


###########reading a file which contains blast hits.
$filenamex = $ARGV[0];
$fileflkfasta = $ARGV[1];
$fileoutputflkseq = $ARGV[2];
$fileoutputlist = $ARGV[3];


#open blast file
unless ( open (PRO, $filenamex) ) {
    print "Cannot open blast result file \"$filenamex\"\n";
    exit;
}
@allalignments = <PRO>;
close PRO;




#find the unique flanking seuqence insertion locations
my %uniquflank; #only record uniqu ID
for my $thisalignment (@allalignments) {
    #look for read information and sequence
    @thisline = split('\s+',$thisalignment); #read ID
    $flankingid=$thisline[0];
    $chr=$thisline[1];
    $alignmentlen=$thisline[2];
    $flankinglen=$thisline[3];
    $chrlen=$thisline[4];
    $flkstart=$thisline[5];
    $flkend=$thisline[6];
    $chrstart=$thisline[7];
    $chrend=$thisline[8];
    $evalue=$thisline[9];
    $bitscore=$thisline[10];

    #determin the insertion location.
    @flkinfo = split('\|',$flankingid); # extract read information
    $readid=$flkinfo[0];
    $targetend=$flkinfo[1];
    $readend=$flkinfo[2];
    @flends=split('-', $flkinfo[4]);
    $smallend=$flends[0];
    $largeend=$flends[1];
    $flkinfo[5]=~m/readlen=(\w*)/;
    $readlength=$1;

    #discard an alignemnt if evalue is too large or if mapped target is too short.
    if( $alignmentlen < $alignment_length_cutoff || $evalue > $evalue_cutoff) {
        #print $thisalignment, "\n";
        next;
    }


    #skip only a small part of flanking seq  was aligned
    if( $largeend-$smallend > $alignmentlen + $diff_alignemnt_total_length ) {
       #  print $thisalignment;
       next;
    }

    #if the mapped part is to far to target end, skip
    if( $targetend eq "target=5" ){
        if ($readend eq "read=5") {
            next if($largeend -$flkend > $gap_to_target_cutoff);
        } else {  #read is 3' end
            next if($flkstart> $gap_to_target_cutoff);
        } 
    } elsif ($targetend eq "target=3"){
        if ($readend eq "read=5") {
            next if($largeend -$flkend> $gap_to_target_cutoff);
        } else { #read is 3' end
            next if($flkstart> $gap_to_target_cutoff);
        } 
    } else {
      print "something must be wrong\n";
      print "$thisalignment\n";
    }

    #determine if it is the unique alignment. 
    #multiple one need to keep the signigicanly best one.
    #do we need do this? --if multiple alignment are the same level good, then keep all alignments.
    if ($uniquealg ne 1){ #not unique alignment keep all alignments
        push (@{$uniquflank{$flankingid}{"alignment"}}, $thisalignment);
    } else {
        if(!exists($uniquflank{$flankingid}) ){
           push (@{$uniquflank{$flankingid}{"alignment"}}, $thisalignment);
           $uniquflank{$flankingid}{"best"}=$thisalignment;
        } else {
           my @oldaln=split('\s+',$uniquflank{$flankingid}{"best"});
           #compare the bitscore with the largest one
           if ($bitscore > $oldaln[10] ){
               push (@{$uniquflank{$flankingid}{"alignment"}}, $thisalignment);
               $uniquflank{$flankingid}{"best"}=$thisalignment;
               #go through all existing alignments and compare with the new best to remove bad ones.
               for ( my $index = scalar @{$uniquflank{$flankingid}{"alignment"}}; $index >= 0; --$index ) {
                     my @oldalnscore=split('\s+', $uniquflank{$flankingid}{"alignment"}[$index]);
                     if ($oldalnscore[10]/$bitscore < $bitscoreratiocutoff ){
                        splice @{$uniquflank{$flankingid}{"alignment"}}, $index, 1;
                     }
               }
               
           } elsif ($bitscore/$oldaln[10] >= $bitscoreratiocutoff ){ 
               #print $uniquread{$readid}{"alignment"};
               #print "-------\n";
               #print $thisalignment;
               #print $uniquflank{$flankingid}{"best"};
               push (@{$uniquflank{$flankingid}{"alignment"}}, $thisalignment);
           }
       }
    }
}

#summarize unique reads used to identification of insertions
my $totallen=0;
my $totalcnt=0;
my %reads=();
foreach my $flkid (sort keys %uniquflank) {
    my @alignments = @{$uniquflank{$flkid}{"alignment"}};
    foreach my $thisalignment (@alignments){
       #print "$thisalignment\n";
       #next;
       @thisline = split('\s+',$thisalignment); #read ID
       $flankingid=$thisline[0];
       @flkinfo = split('\|',$flankingid); # extract read information
       $readid=$flkinfo[0];
       $flkinfo[5]=~m/readlen=(\w*)/;
       $readlength=$1;
       
       if(!exists($reads{$readid})){
          $totallen=$totallen+$readlength;
          $totalcnt++;
          $reads{$readid}=$readid;
       }
    }
}

$avglen=int($totallen/$totalcnt);
print "Total number of used reads:$totalcnt\n";
print "Average length of used reads:$avglen\n";


### count all flanking seqs for each insertion location
my %inserts; #unqi insertion location
foreach my $flkid (sort keys %uniquflank) {
    #look for read information and sequence
    #if(!exists($uniquflank{$flkid}{"alignment"})){
    #   print "$flkid\n";
    #}
    # my $thisalignment = $uniquflank{$flkid}{"alignment"};
    my @alignments = @{$uniquflank{$flkid}{"alignment"}};
    foreach my $thisalignment (@alignments){
       #print "$thisalignment\n";
       #next;
       @thisline = split('\s+',$thisalignment); #read ID
       $flankingid=$thisline[0];
       $chr=$thisline[1];
       $alignmentlen=$thisline[2];
       $flankinglen=$thisline[3];
       $chrlen=$thisline[4];
       $flkstart=$thisline[5];
       $flkend=$thisline[6];
       $chrstart=$thisline[7];
       $chrend=$thisline[8];
       $evalue=$thisline[9];
       $bitscore=$thisline[10];
       
       #print $thisalignment;
       #next;
       
       #determin the insertion location.
       @flkinfo = split('\|',$flankingid); # extract read information
       $readid=$flkinfo[0];
       $targetend=$flkinfo[1];
       $readend=$flkinfo[2];
       @flends=split('-', $flkinfo[4]);
       $smallend=$flends[0];
       $largeend=$flends[1];
       
       #if the flanking sequence is on the 3' of target sequence,
       #then if (flanking seqence is on 3' of read) , then the small end is insertion .
       #     if (flanking seqence is on 5' of read) , then the large end is insertion 
       my $location;
       if($readend eq "read=3"){
          $location=$chrstart;
       } elsif ($readend eq "read=5"){
          $location= $chrend;
       } else {
         print "something must be wrong\n";
       }
      
       $locationbin=int($location/$bin_size);
       if(exists($inserts{$chr}{$locationbin})){
          if($targetend eq "target=5"){
             if (!exists($inserts{$chr}{$locationbin}{"count5"})){
                 $inserts{$chr}{$locationbin}{"count5"}=1;
                 $inserts{$chr}{$locationbin}{"cut5"}=$location;
                 $inserts{$chr}{$locationbin}{"cut5read"}=$thisalignment;
                 push(@{$inserts{$chr}{$locationbin}{"alignments5"}}, $thisalignemnt);
             } else {
                 $inserts{$chr}{$locationbin}{"count5"}=$inserts{$chr}{$locationbin}{"count5"}+1;
                 push (@{$inserts{$chr}{$locationbin}{"alignments5"}}, $thisalignemnt);
                 if($readend eq "read=3"){
                    if($chrstart < $chrend && $inserts{$chr}{$locationbin}{"cut5"} > $chrstart){
                       $inserts{$chr}{$locationbin}{"cut5"}=$chrstart;
                       $inserts{$chr}{$locationbin}{"cut5read"}=$thisalignment;
                    }
                    if($chrstart > $chrend && $inserts{$chr}{$locationbin}{"cut5"} < $chrstart){
                       $inserts{$chr}{$locationbin}{"cut5"}=$chrstart;
                       $inserts{$chr}{$locationbin}{"cut5read"}=$thisalignment;
                    }
                 } elsif ($readend eq "read=5"){
                    $location=$chrend;
                    if($chrstart < $chrend && $inserts{$chr}{$locationbin}{"cut5"} < $chrend){
                       $inserts{$chr}{$locationbin}{"cut5"}=$chrend;
                       $inserts{$chr}{$locationbin}{"cut5read"}=$thisalignment;
                    }
                    if($chrstart > $chrend && $inserts{$chr}{$locationbin}{"cut5"} > $chrend){
                       $inserts{$chr}{$locationbin}{"cut5"}=$chrend;
                       $inserts{$chr}{$locationbin}{"cut5read"}=$thisalignment;
                    }
                 } else {
                    print "something must be wrong\n";
                 }
             }
          } elsif ($targetend eq "target=3") {
             if (!exists($inserts{$chr}{$locationbin}{"count3"})){
                 $inserts{$chr}{$locationbin}{"count3"}=1;
                 $inserts{$chr}{$locationbin}{"cut3"}=$location;
                 $inserts{$chr}{$locationbin}{"cut3read"}=$thisalignment;
                 push(@{$inserts{$chr}{$locationbin}{"alignments3"}}, $thisalignemnt);
             } else {
                 $inserts{$chr}{$locationbin}{"count3"}=$inserts{$chr}{$locationbin}{"count3"}+1;
                 push(@{$inserts{$chr}{$locationbin}{"alignments3"}}, $thisalignemnt);
                 if($readend eq "read=3"){
                    if($chrstart < $chrend && $inserts{$chr}{$locationbin}{"cut3"} > $chrstart){
                       $inserts{$chr}{$locationbin}{"cut3"}=$chrstart;
                       $inserts{$chr}{$locationbin}{"cut3read"}=$thisalignment;
                    }
                    if($chrstart > $chrend && $inserts{$chr}{$locationbin}{"cut3"} < $chrstart){
                       $inserts{$chr}{$locationbin}{"cut3"}=$chrstart;
                       $inserts{$chr}{$locationbin}{"cut3read"}=$thisalignment;
                    }
                 } elsif ($readend eq "read=5"){
                    $location=$chrend;
                    if($chrstart < $chrend && $inserts{$chr}{$locationbin}{"cut3"} < $chrend){
                       $inserts{$chr}{$locationbin}{"cut3"}=$chrend;
                       $inserts{$chr}{$locationbin}{"cut3read"}=$thisalignment;
                    }
                    if($chrstart > $chrend && $inserts{$chr}{$locationbin}{"cut3"} > $chrend){
                       $inserts{$chr}{$locationbin}{"cut3"}=$chrend;
                       $inserts{$chr}{$locationbin}{"cut3read"}=$thisalignment;
                    }
                 } else {
                    print "something must be wrong\n";
                 }
             }
          } else {
             print "something must be wrong\n";
          }
       } else {
          if ($targetend eq "target=3"){
              $inserts{$chr}{$locationbin}{"count3"}=1;
              $inserts{$chr}{$locationbin}{"cut3"}=$location;
              $inserts{$chr}{$locationbin}{"cut3read"}=$thisalignment;
              push(@{$inserts{$chr}{$locationbin}{"alignments3"}}, $thisalignemnt);
          }
          if ($targetend eq "target=5"){
              $inserts{$chr}{$locationbin}{"count5"}=1;
              $inserts{$chr}{$locationbin}{"cut5"}=$location;
              $inserts{$chr}{$locationbin}{"cut5read"}=$thisalignment;
              push(@{$inserts{$chr}{$locationbin}{"alignments5"}}, $thisalignemnt);
          }
       }
   }
}


###########reading gtf file.
#$filenames = "ref.gtf";
$gtffilename = $ENV{'gtffilename'};
unless ( open (PRO, $gtffilename) ) {print "Cannot open file \"$gtffilename\"\n\n";exit;}
@gtf = <PRO>;
close PRO;
my %allgens=();
for ($g=0; $g<= $#gtf; $g++) {
    @fs = split ("\t", $gtf[$g]);
    $chr=$fs[0];
    $gene=$fs[2];
    next if ($gene ne "gene");
    $gstart=$fs[3];
    $gend=$fs[4];
    $genedescript=$fs[8];
    if (!$genedescript){print "$gtf[$g]\n";next;}
    #$genedescript =~ m/gene_id\s"(\w*)";\s/;
    #$genedescript =~ m/gene_id\s"(\w*)";\s/; #for GTF
    #$genedescript =~ /gene_id "(.*)\"\; /; #for GFF3
    @spg = split('gene_id \"',$genedescript);
    @spg2 = split('\"',$spg[1]);
    $geneid = $spg2[0];
    $posbin=int($gstart/$genome_bin);
    $allgene{$chr}{$posbin}{$geneid}{"genestart"}=$gstart;
    $allgene{$chr}{$posbin}{$geneid}{"geneend"}=$gend;
    #CZ 3/19/2021
    if( int($gstart/$genome_bin) ne int($gend/$genome_bin) ){
        $posbin=int($gend/$genome_bin);
        $allgene{$chr}{$posbin}{$geneid}{"genestart"}=$gstart;
        $allgene{$chr}{$posbin}{$geneid}{"geneend"}=$gend;
    }
#print "$posbin\t $geneid \t $chr \t $allgene{$chr}{$posbin}{$geneid}{'genestart'} \t $allgene{$chr}{$posbin}{$geneid}{'geneend'}\n";
}


#open falnking seq fasta file
unless ( open (PRO, $fileflkfasta) ) {
    print "Cannot open read fasta file \"$ARGV[1]\"\n";
    exit;
}
@allfaseq = <PRO>;
close PRO;
my %allflkseq=();
for ($f=0; $f<= $#allfaseq; $f++) {
    if ($allfaseq[$f] =~ m/>/) {
        chomp $allfaseq[$f];
        @forindiseq = split('>',$allfaseq[$f]);
        @forindiseq2 = split(' ',$forindiseq[1]);
        $nextline = $f+1;
        $nextseq = $allfaseq[$nextline];
        chomp $nextseq;
        $allflkseq{$forindiseq2[0]}=$nextseq;
    }
}


#open output falnking seq 
unless ( open (OUTFLKSEQ, '>', $fileoutputflkseq) ) {
    print "Cannot open to write file \"$fileoutputflkseq\"\n";
    exit;
}

#open output insertion list file
unless ( open (OUTLIST, '>', $fileoutputlist) ) {
    print "Cannot open to write file \"$fileoutputlist\"\n";
    exit;
}

#display all inserts 
foreach my $chrr (sort keys %inserts) {
   foreach my $cord (sort keys %{$inserts{$chrr}}) {
       if(!exists( $inserts{$chrr}{$cord}{"cut3"})){ 
          $ct3="XXX";
          $cnt3=0;
          if(!exists( $inserts{$chrr}{$cord}{"cut5"})){
             print "something must be wrong\n";
          } else {
             $ct5=$inserts{$chrr}{$cord}{"cut5"};
             $cnt5=$inserts{$chrr}{$cord}{"count5"};
             $insertstart=$ct5;
             $insertend = $ct5;
          }
       } else {
          $ct3=$inserts{$chrr}{$cord}{"cut3"};
          $cnt3=$inserts{$chrr}{$cord}{"count3"};
          if(!exists( $inserts{$chrr}{$cord}{"cut5"})){
            $ct5="XXX";
            $cnt5=0;
            $insertstart=$ct3;
            $insertend = $ct3;
          } else {
             $ct5=$inserts{$chrr}{$cord}{"cut5"};
             $cnt5=$inserts{$chrr}{$cord}{"count5"};
             if($ct3 >= $ct5){
                $insertstart= $ct5;
          	$insertend = $ct3;
       	     } else {
          	$insertstart= $ct3;
                $insertend = $ct5;
             }  
          }
       }

       my $realloc=$insertstart."-".$insertend;
       my $insertgap=$insertend-$insertstart;
       my $totalcount=$cnt3+$cnt5;



       $genebin=int($insertend/$genome_bin);
       $findone=0;
       foreach my $geneid (sort keys %{$allgene{$chrr}{$genebin}})
       {
#print "Chr$chrr:$realloc\t $geneid\t $allgene{$chrr}{$genebin}{$geneid}{'genestart'} \t $allgene{$chrr}{$genebin}{$geneid}{'geneend'}\n";
           if( $allgene{$chrr}{$genebin}{$geneid}{"genestart"} - 2000 <=$insertstart &&
               $allgene{$chrr}{$genebin}{$geneid}{"geneend"} + 2000 >= $insertend )
           {
	       $geneinf=$geneid."|"."Chr".$chrr.":".$allgene{$chrr}{$genebin}{$geneid}{"genestart"}."-".$allgene{$chrr}{$genebin}{$geneid}{"geneend"};
               $findone=1;
               goto mybreak;
           } 
       }
       mybreak:
       if($findone eq 0)
       {
          $geneinf = "Intergenic";
       }
       print OUTLIST "Chr$chrr:$realloc\t $insertgap\t $totalcount\t $geneinf\n";

       #output to a file for flanking sequences
       print OUTFLKSEQ "\n\n****Chr$chrr:$realloc\t $insertgap\t $totalcount\t $geneinf\n";
       if(exists($inserts{$chrr}{$cord}{"cut5read"})){
          @thisline = split('\s+',$inserts{$chrr}{$cord}{"cut5read"}); #read ID
          $flankingid=$thisline[0];
          print OUTFLKSEQ ">5' of the target|$inserts{$chrr}{$cord}{'cut5read'}";

          #determin the read direction. get complement reverse seq if need
          @flkinfo = split('\|',$flankingid); # extract read information
          $flkinfo[1]=~m/target=(\w*)/;
          $targetend=$1;
          $flkinfo[2]=~m/read=(\w*)/;
          $readend=$1;
          if($targetend ne $flkinfo[2]){
             my $revcomp = reverse $allflkseq{$flankingid};
             $revcomp =~ tr/ATGCatgc/TACGtacg/;
             $allflkseq{$flankingid} = $revcomp ;
          }
          print OUTFLKSEQ "$allflkseq{$flankingid}\n";
       }
       if(exists($inserts{$chrr}{$cord}{"cut3read"})){
          @thisline = split('\s+',$inserts{$chrr}{$cord}{"cut3read"}); #read ID
          $flankingid=$thisline[0];
          print OUTFLKSEQ ">3' of the target|$inserts{$chrr}{$cord}{'cut3read'}";

          #determin the read direction. get complement reverse seq if need
          @flkinfo = split('\|',$flankingid); # extract read information
          $flkinfo[1]=~m/target=(\w*)/;
          $targetend=$1;
          $flkinfo[2]=~m/read=(\w*)/;
          $readend=$1;
          if($targetend ne $flkinfo[2]){
             my $revcomp = reverse $allflkseq{$flankingid};
             $revcomp =~ tr/ATGCatgc/TACGtacg/;
             $allflkseq{$flankingid} = $revcomp ;
          }
          print OUTFLKSEQ "$allflkseq{$flankingid}\n";
       }

   }
}
close(OUTFLKSEQ);
close(OUTLIST);
