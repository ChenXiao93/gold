#!/bin/bash
for var in    "400" "500" "700" 
do
    rm ./$var/inp
    rm ./$var/run
    #rm ./$var/slurm*

    cp ./inp.kick  ./$var/inp
    cp ./run  ./$var/run

    sed -i "s/delay/$var/g"  ./$var/inp
    sed -i "s/delay/$var/g"  ./$var/run
    
    cd ./$var/
    sbatch run
    cd ..
done
