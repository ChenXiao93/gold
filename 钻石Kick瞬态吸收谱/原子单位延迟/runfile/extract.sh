#!/bin/bash
#for var in "108" "477"  "518" "559" "600"  "641"  "682" "723" "1092"  "pump" "probe"
for var in "100" "200" "300" "400" "500" "600" "700" "800" "900" "1000" "1100" "1200"

do
    cp -r ./$var/td.general  ./newdata/$var
done
