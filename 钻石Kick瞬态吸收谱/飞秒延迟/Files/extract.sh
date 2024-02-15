#!/bin/bash
for var in "108" "477"  "518" "559" "600"  "641"  "682" "723" "1092"  "pump" "probe"
do
    cp -r ./$var/td.general  ./data/$var
done
