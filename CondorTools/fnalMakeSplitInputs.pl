#!/usr/bin/perl

if ($#ARGV < 0){
   print "\nUsage: splitInput numFilePerSnippet\n";
   print "split the inputList.txt into files with numFilePerSnippet each. \n\n";
   exit;
}

print "$ARGV[0]  $ARGV[1]\n";
$numArgs = $#ARGV + 1;
print "number of args = $numArgs\n";

#if ($#ARGV == 1){
if ($numArgs == 2){
   
   $perNumFile = $ARGV[0];
   $nJobs = $ARGV[1];

   print "nJobs      : $nJobs\n";
   print "perNumFile : $perNumFile\n";

   #open(IN, "wjets_wz.list");
   #open(IN, "qcd_pt1800.list");
   open(IN, "file.list");
   $cnt = -1;
   $appendix = 0;

   while(<IN>){
      $cnt ++;
      if( $cnt % $perNumFile ==0 ){
         if( $cnt !=0 ){ close(OUT); }
#         $outFileName = sprintf("inputList_%03d.txt", $appendix);
         $outFileName = sprintf("inputList_%d.rtlst", $appendix);
         open(OUT, ">$outFileName");
         $appendix++;
      }
      print OUT $_;
		last if $appendix > $nJobs
   }
   close(OUT);
}

