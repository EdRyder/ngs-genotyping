#!/usr/bin/perl

use Bio::Seq;
use Bio::Seq::Quality;
use Bio::SeqIO;
use File::Slurp;
use DDP;
use DBI;
use List::MoreUtils qw(uniq);
use strict;

###
#   USAGE: perl crispr_miseq-scanner-fastq-batch.pl <results-file> <debug-file> <run_id> 
#
###

## CHANGE THESE SETTINGS FOR YOUR SERVER
my $dsn = "DBI:mysql:mouseNGS";
my $username = "root"; ## 
my $password = "";
#my $username = "nobody";
#my $password = "";

my $dbh = DBI->connect($dsn,$username,$password);


## GET FILES AND PUT THEM INTO AN ARRAY, THEN USE THAT ARRAY TO GET R1 and R2 paired ends for each sample. I'm sure there's an easier way of doing this...
my @name_array = ();
my $dir = ".";
my @files = grep ( -f ,<*.fa>);
foreach my $files(@files){
   my @samples = split (/_/, $files);
   my $name = $samples[0];
   push @name_array, $name;
}
my @uniq_name = uniq(@name_array);
##

open (OUT, ">$ARGV[1]");
print (OUT "Bin	Sequence	score\n");

my $filetemp = $ARGV[0];
open (OUT2, ">", $filetemp);
print (OUT2 "Mouse	Genotype	Gene	Plate	Comments\n");

#$run_id = "run133";
my $run_id = $ARGV[2];

## get input file of crisprs and read into an array
##@crisprs = read_file($crisprs_in, chomp => 1);

foreach my $uniq_name(@uniq_name){
   my @crisprs = ();
   my $genotype = "";
   my @wt_count_array = ();
   my @mut_count_array = ();
   my  $yy = "";
   my   $mouse = "" ;
   my   $plate = "" ;
   my  $gene = "" ;

   #GET CRISPR INFO FROM DATABASE
   my $sth = $dbh->prepare("select distinct t1.*, t2.*  from miseq_mice t1 LEFT JOIN miseq_sequences t2 on t1.gene_id = t2.gene_id where t1.seq_id = '$uniq_name' and plate_id = '$run_id' and match_strand = 'plus' ");
   $sth->execute();
   while (my @row = $sth->fetchrow_array){
     #print "$row[1], $row[2], $row[3], $row[4], $row[7], $row[9]\n";
        $yy = "$row[7]\t$row[9]";
        $mouse = $row[2];
        $plate = $row[3];
        $gene = $row[4];
     push @crisprs, $yy;
   }
   #
 
   my $grepmatch = "$uniq_name"."_\*_R1_\*\.fa";
   my @files2 = grep(-f ,<$uniq_name\_*_R1_*.fa>); #only looking at one end for the test. Remove R1 to get both ends.
   #my @files2 = grep ( -f ,<$grepmatch>); #only looking at one end for the test. Remove R1 to get both ends.
   ## Process each file and look for matches in the sequence to the list of crisprs 
   foreach my $file (@files2) {
      my $r = 0;
      my $tt = 25; #threshold score
      print "processing $file uniq_name $uniq_name for  mouse $mouse $plate  $gene grepmatch = $grepmatch  \n";
      #p(@crisprs);
      my $testseq = Bio::SeqIO->new ('-format' => 'fasta',
                            '-file' => $file
                           );
      my @mut_array = ();
      my @wt_array = ();
      my @other_array = ();
      my @filter_array = ();
      my @contamination_array = ();
      while (my $s = $testseq->next_seq() ){ ## add in loop to check all crisprs in input file
         #$s->threshold($tt); 
         #p($s);
         my  $count = 149 ;# left over from other version
         #check quality
         #my $quals = $s->qual_text;
         #print "quals = $quals\n";
         #$count = $s->clear_ranges_length($tt);
         if ($count > 145){ # only include sequences of high quality throughout 150bp run
            #print "$count\n";
            #

           my $sequence = $s->seq();
           my $hit = 0;
           my $cont = 0;
           #print ( "sequence is ", $sequence, "\n");
           if ($sequence =~ /GAGCTGTACAAGTAAGGTGTGGGAGGTTTTAGATCGGAAGAGCGGTTCAGCA/){ 
              push (@contamination_array, $sequence);
          #    print OUT "contamination	$sequence	$count\n";
              $cont = 1;
           }
           else{
             foreach my $crispr (@crisprs){
                my @crispr_split = split(/\t/, $crispr);
                my $crispr_off = uc($crispr_split[1]);
                my $crispr_off_lc = $crispr_split[1];
                my $crispr_rev = reverse($crispr_off);
                my $crispr_rev =~ tr/ACGTacgt/TGCAtgca/;
                my $bin = $crispr_split[0];
                #print "sequence = $sequence, crispr_off = $crispr_off, crispr_rev  = $crispr_rev, bin = $bin, hit = $hit\n";
                if ($sequence =~ /$crispr_off/ && $hit == 0){
                   if ($bin eq "filter"){
                      push (@filter_array, $sequence);
           #           print OUT "filter	$sequence	$count	$r\n";
                      $hit++;
                   }
                   elsif ($bin eq "mut"){
                      push (@mut_array, $sequence);
           #           print OUT "mut	$sequence	$count	$r\n";
                      $hit++;
                   }
                   elsif ($bin eq "wt"){
                      push (@wt_array, $sequence);
           #           print OUT "wt	$sequence	$count	$r\n";
                      $hit++;
                   }
                }
             }
             if ($hit == 0 && $cont == 0){
                push (@other_array, $sequence);
          #      print OUT "other	$sequence	$count\n";
             }
          }
          $r++;
        }
     }       
     my $count_mut = @mut_array;
     push @mut_count_array, $count_mut;
     my $count_wt = @wt_array;
     push @wt_count_array, $count_wt;
     print "mut array = ".@mut_array."\n";
     print "wt array = ".@wt_array."\n";
     print "other array = ".@other_array."\n";
     print "contamination_array = ".@contamination_array."\n";
     print "filter array = ".@filter_array."\n";
     my $filtered_wt = @wt_array - @filter_array;
     print "wt - filtered = $filtered_wt\n";
   } # end of each file
   
   ## GENOTYPING PROCESSING INFORMATION
   #@wt_count_array = ();
   #@mut_count_array = ();

   my $avg_mut = average(@mut_count_array);
   my $avg_wt = average(@wt_count_array);
   if ($avg_wt + $avg_mut == 0){
      $avg_wt = 1;
   }
   print "avg_mut = $avg_mut, avg_wt = $avg_wt\n";
   my $ratio = sprintf("%.2f", ($avg_wt/($avg_wt + $avg_mut)));
   print "ratio = $ratio\n";
   if ($ratio < 0.15 && $avg_mut > 1000){
      $genotype = "Hom";
   }
   elsif ($ratio < 0.75 && $ratio > 0.25 && $avg_mut > 1000 && $avg_wt > 1000){
      $genotype = "Het";
   }
   elsif ($ratio > 0.85 && $avg_wt > 1000){
       $genotype = "WT";
   } 
   else {
      $genotype = "retest";
   }
   my $comments = "WT count = $avg_wt, Mut count = $avg_mut, Ratio = $ratio";
   print (OUT2 "$mouse	$genotype	$gene	$plate	$comments\n");
   ##
} # end of all files
close OUT;
close OUT2;


sub average {
   my @array = @_; # save the array passed to this function
   my $sum; # create a variable to hold the sum of the array's values
   foreach (@array) { $sum += $_; } # add each element of the array 
   # to the sum
   return $sum/@array; # divide sum by the number of elements in the
   # array to find the mean eg
   # @dataArray = (1, 2, 3, 4, 5, 6, 7, 8, 9);
   # $avgOfArray = average(@dataArray);
   # here the average will be 5

}

