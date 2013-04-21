#!/usr/bin/perl

if ($#ARGV < 1){
   print "\nUsage: $0 numOfFilePerJob fileList.txt\n";
   print "split the fileList.txt into files with numOfFilesPerJob each. \n\n";
   exit;
}

if($#ARGV ==1 ){
   
   $iFilesPerJob = $ARGV[0];
   $sFileList = $ARGV[1];
   print "iFilesPerJob : $iFilesPerJob\n";
   print "sFileList : $sFileList\n";

	#use the same file name for input lists
	#if there is no '.' seperated extension, use the whole name
	#$iDotIndex = rindex($sFileList)
	#if ($iDotIndex > -1) {
	#	$sName = subtr($sFileList,0, $iDotIndex);
	#} else {
	#	$sName = $sFileList;
	#}
	print "Creating input file list with name $sName";


   #open(IN, "file.list");
   open(IN, $sFileList);
   $cnt = -1;
   $appendix = 0;

   while(<IN>){
      $cnt ++;
      if( $cnt % $iFilesPerJob ==0 ){
         if( $cnt !=0 ){ close(OUT); }
#         $outFileName = sprintf("inputList_%03d.txt", $appendix);
         $outFileName = sprintf("inputList_%d.rtlst", $appendix);
         #$outFileName = sprintf("%s_%d.rtlst",$sName, $appendix);
         open(OUT, ">$outFileName");
         $appendix++;
      }
      print OUT $_;
   }
   close(OUT);
}

