#!/bin/csh
#release jobs on holds on condor

set firstjob = 17181
set lastjob = 17358 

while ($firstjob <= $lastjob)
	condor_release $firstjob
	@ firstjob++
end

