#!/bin/csh
#holds jobs running on condor

set firstjob = 74432
set lastjob = 74483

while ($firstjob <= $lastjob)
	condor_rm $firstjob
	@ firstjob++
end

