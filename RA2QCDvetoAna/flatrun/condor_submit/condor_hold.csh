#!/bin/csh
#holds jobs running on condor

set firstjob = 17180
set lastjob = 17358 

while ($firstjob <= $lastjob)
	@ firstjob++
	condor_hold $firstjob
end

