#! /bin/csh

#makes files list used to run on FIU farm and local
set srcDir = "/share/store/users/samantha/11282012_Ntuples_chsJetswithFilter"
set desDir = "/mnt/hadoop/cms/store/user/samantha/11282012_Ntuples_chsJetswithFilter"

set dirs = "QCD1 QCD2 QCD3 QCD4 QCD5 QCD6 QCD7"

foreach dir ( $dirs )

	set compdir = "$srcDir/$dir"
	if ( -e $compdir ) then
		echo "Processing $compdir"
		set files = `ls $compdir`
		set outfile = "$srcDir/${dir}_wFilters_farm.files"
		set outfile_local = "$srcDir/${dir}_wFilters_local.files"
		echo $files
		echo $outfile
		echo $outfile_local
		if ( -e $outfile ) then
			echo "Output file with the name $outfile already exists! plean remove it. otherwise it may be appended with duplicate files!!!!"
			exit
		endif
		if ( -e $outfile_local ) then
			echo "Output file with the name $outfile_local already exists! plean remove it. otherwise it may be appended with duplicate files!!!!"
			exit
		endif

		touch $outfile
		touch $outfile_local
		foreach file ( $files )
			echo "writin file $file"
			echo "${desDir}/${dir}/$file" >> $outfile 
			echo "${srcDir}/${dir}/$file" >> $outfile_local 
		end
	else 
		echo "$dir does not exist. skipping.."
	endif
end
