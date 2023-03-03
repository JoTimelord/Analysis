#!/bin/bash

declare -a FILENAMES=("2oslep_2ak4_1ak8_v2"
                      "2oslep_2ak4_1ak8" 
                      "1lep_1ak8_2ak4_v1")

rm -r ~/Analysis/outputs
mkdir ~/Analysis/outputs

for i in "${FILENAMES[@]}"; do
    echo "$i"
    rm -r /home/users/joytzphysics/Analysis/outputs/output_semiMerge_"$i"
    rm -r /home/users/joytzphysics/Analysis/outputs/output_fullyMerge_"$i"
    rm -r /home/users/joytzphysics/Analysis/outputs/semiMerge_logfiles_"$i"
    rm -r /home/users/joytzphysics/Analysis/outputs/fullyMerge_logfiles_"$i"
    mkdir /home/users/joytzphysics/Analysis/outputs/output_semiMerge_"$i"
    mkdir /home/users/joytzphysics/Analysis/outputs/output_fullyMerge_"$i"
    mkdir /home/users/joytzphysics/Analysis/outputs/semiMerge_logfiles_"$i"
    mkdir /home/users/joytzphysics/Analysis/outputs/fullyMerge_logfiles_"$i"
done




