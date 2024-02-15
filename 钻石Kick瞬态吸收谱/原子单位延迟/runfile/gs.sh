#!/bin/bash
for var in "100" "200"  "300" "400" "500"  "700"  "800" "900" "1000" "1100" "1200"
do
    mkdir $var
    
    cp ./inp.gs  ./$var/inp
    cp ./run.fast  ./$var/run

    cd ./$var/
    sbatch run
    cd ..
done
