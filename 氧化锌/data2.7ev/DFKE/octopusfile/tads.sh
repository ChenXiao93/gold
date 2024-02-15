#!/bin/bash
#for var in  "17" "30" "43"  "56" "97" # "137" "177"
for var in "137" "177"
do
     #mkdir  A$var
     cp ./inp.kick  ./A$var/inp
     cp ./run	./A$var/run

     sed -i "s/amp/$var/g"  ./A$var/inp
     sed -i "s/amp/$var/g"  ./A$var/run

done
