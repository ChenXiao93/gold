#!/bin/bash
for var in "108" "477"  "518" "559" "600"  "641"  "682" "723" "1092"  "pump" "probe"
do
#    mkdir $var
    rm ./$var/slurm*
    rm ./$var/inp
    rm ./$var/run
    rm -r ./$var/td.general
    rm -r ./$var/static
    
    cp ./inp.gs  ./$var/inp
    cp ./run.fast  ./$var/run

    cd ./$var/
    sbatch run
    cd ..
done
