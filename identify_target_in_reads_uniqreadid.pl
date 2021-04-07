#!/usr/bin/perl

# input one blast file against target.fa, and further blast against soybean genome for coordinates and inserted genes.

use warnings;


#mapped target length cutoff
#$targetmappingcutoff=30;
$targetmappingcutoff=$ENV{'targetmappingcutoff'};

#alignment evalu cutoff
#$evaluecutoff=0.0000000001;
$evaluecutoff=$ENV{'evaluecutoff'};

#alignment of target missed end length tolerence. if one end missing too long seq, don't consider this end of an alignment .
#$target_end_error_tolerence = 9;
$target_end_error_tolerence = $ENV{'target_end_error_tolerence'};

#flanking seuqnece length cutoff
#$flanking_length_cutoff = 50;    
$flanking_length_cutoff = $ENV{'flanking_length_cutoff'};    

###########reading a file which contains blast hits.
$filenamex = $ARGV[0];


#open read fasta file
unless ( open (PRO, $ARGV[1]) ) {
    print "Cannot open read fasta file \"$ARGV[1]\"\n";
    exit;
}
@allfaseq = <PRO>;
close PRO;
for ($f=0; $f<= $#allfaseq; $f++) {
    if ($allfaseq[$f] =~ m/>/) {
        chomp $allfaseq[$f];
        @forindiseq = split('>',$allfaseq[$f]);
        @forindiseq2 = split(' ',$forindiseq[1]);
        $nextline = $f+1;
        $nextseq = $allfaseq[$nextline];
        chomp $nextseq;
        $allfaseq{$forindiseq2[0]}=$nextseq;
    }
}

#open blast file
my @allalignemnts = ();
unless ( open (PRO, $filenamex) ) {
    print "Cannot open blast result file \"$filenamex\"\n";
    exit;
}
@allalignments = <PRO>;
close PRO;

my %uniqueread; #only record uniqu ID
for $thisalignment (@allalignments) {

    #look for read information and sequence
    @thisline = split('\s+',$thisalignment); #read ID
    $readid=$thisline[1];

    $targetlen=$thisline[3];
    $readlen=$thisline[4];
    $targetstart=$thisline[5];
    $targetend=$thisline[6];
    $readstart=$thisline[7];
    $readend=$thisline[8];
    $evalue=$thisline[9];
    $bitscore=$thisline[10];

 
    #don't need to identify the unique id.
    if(exists($uniquread{$readid})){
       my @oldaln=split('\s+',$uniquread{$readid}{"alignment"});
       if ($oldaln[10] > $bitscore ){ #keep the one that has the large bitscore
          #print $uniquread{$readid}{"alignment"};
          #print $thisalignment;
          #next;
       }
    };

    #discard an alignemnt if evalue is too large or if mapped target is too short.
    if( $targetend-$targetstart < $targetmappingcutoff || $evalue > $evaluecutoff) {
        # print $thisalignment, "\n";
        next;
    }

    #discard an alignemnt if target both ends could NOT  be mapped acurately.
    if($targetlen-$targetend >= $target_end_error_tolerence && $targetstart >= $target_end_error_tolerence ) {
       #print $thisalignment;
       next;
    }

    #print $thisalignment;
    #next;

    #determin the 3' flanking sequence .
    if($targetlen-$targetend < $target_end_error_tolerence ) {
       #print $thisalignment;
       # print "$readid \t $readlen \t $readstart \t $readend\n";

       #determine if the mapping direction of reads
       my $title; my $flankseq;
       if ($readend > $readstart){
           $flength=$readlen-$readend;
           next if($flength < $flanking_length_cutoff); #discard if flanking sequence is too short
           $title= ">$readid|target=3|read=3|length=$flength|$readend-$readlen|readlen=$readlen|target=$targetlen-$targetstart-$targetend";       
           $flankseq=substr($allfaseq{$readid}, $readend);
           #print "$flankseq\n";
           $uniquread{$readid}{"alignment"} = $thisalignment;
           $uniquread{$readid}{"title3"} = $title;
           $uniquread{$readid}{"flankseq3"} = $flankseq;
       } else {
           next if($readend < $flanking_length_cutoff); #discard if flanking sequence is too short
           $title= ">$readid|target=3|read=5|length=$readend|1-$readend|readlen=$readlen|target=$targetlen-$targetstart-$targetend";  
           $flankseq=substr($allfaseq{$readid},1, $readend);
           #print "$flankseq\n";
           $uniquread{$readid}{"alignment"} = $thisalignment;
           $uniquread{$readid}{"title3"} = $title;
           $uniquread{$readid}{"flankseq3"} = $flankseq;
       }
    }
    #determin the 5' flanking sequence .
    if($targetstart < $target_end_error_tolerence ) {
       #determine if the mapping direction of reads
       if ($readend > $readstart){
           next if($readstart < $flanking_length_cutoff); #discard if flanking sequence is too short
           $title= ">$readid|target=5|read=5|length=$readstart|1-$readstart|readlen=$readlen|target=$targetlen-$targetstart-$targetend";
           $flankseq=substr($allfaseq{$readid}, 1, $readstart);
           #print "$flankseq\n";
           $uniquread{$readid}{"alignment"} = $thisalignment;
           $uniquread{$readid}{"title5"} = $title;
           $uniquread{$readid}{"flankseq5"} = $flankseq;

       } else {
           $flength=$readlen-$readstart;
           next if($flength < $flanking_length_cutoff); #discard if flanking sequence is too short
           $title=">$readid|target=5|read=3|length=$flength|$readstart-$readlen|readlen=$readlen|target=$targetlen-$targetstart-$targetend";
           $flankseq=substr($allfaseq{$readid}, $readstart);
           #print "$flankseq\n";
           $uniquread{$readid}{"alignment"} = $thisalignment;
           $uniquread{$readid}{"title5"} = $title;
           $uniquread{$readid}{"flankseq5"} = $flankseq;
       }       
    }


}

foreach my $rid (sort keys %uniquread) {
   if(exists($uniquread{$rid}{'title5'})) {
      print "$uniquread{$rid}{'title5'}\n";
      print "$uniquread{$rid}{'flankseq5'}\n";
   }
   if(exists($uniquread{$rid}{'title3'})) {
      print "$uniquread{$rid}{'title3'}\n";
      print "$uniquread{$rid}{'flankseq3'}\n";
   }
}
