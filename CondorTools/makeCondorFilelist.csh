#! /bin/csh

if ( $#argv < 1) then
	echo "Need absolute path as input"
	exit
endif

set dir = $argv[1]
if ( ! -d $dir ) then
	echo "Dir not found! $dir"
	exit
endif

#this is to avoid getting difference anwer if a '/' is put at the end of the path
#set lastpart1 = `echo $dir | awk -F "/" '{print $(NF-1)}'`
set lastpart2 = `echo $dir | awk -F "/" '{print $(NF)}'`
#set lastchr=`${dir#${dir%?}}`
#echo "lastchar = $lastchr "
#	echo $lastpart1
	echo $lastpart2
#if ( $lastpart1 != $lastpart2 ) then
#	echo "lastpart names do not match. Try removing the last '/' charactor!"
#	exit
#else
#	echo "lastpart names match"
#endif

set fileName = $lastpart2".list"
echo "fileName = $fileName"

set files = `ls -1 ${dir}/*.root`

if ( $#files <1 ) then
	echo "No root files found in the dir!"
	exit
endif

if (-e $fileName ) then
	echo "Output file with the name $file already exist. Overwrite it (Y/N)?"
	set input = $<
	if ( $input == "Y" || $input == "y" ) then
		rm -vf $fileName
	else
		exit
	endif
endif

foreach file ( $files )
	echo "dcap:/$file" >> $fileName
end
