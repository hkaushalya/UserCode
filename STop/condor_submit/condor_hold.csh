#!/bin/csh
#holds jobs running on condor

set firstjob = 55255
set lastjob = 56590

while ($firstjob <= $lastjob)
        condor_hold $firstjob
        @ firstjob++
end
