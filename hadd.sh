#!/bin/bash

PROCESSES="DYJetsToLL \
WZH \
TTTo2L2Nu \
WWTo2L2Nu"

rm -r hadds/
mkdir -p hadds/

rm -f .hadd.txt

for PROCESS in ${PROCESSES}; do
    FILENAMES=(${PROCESS}*.root)
    echo "hadd -f hadds/${PROCESS}.root output/${FILENAMES[@]} > hadds/${PROCESS}.log 2>&1" >> .hadd.txt
done

xargs.sh .hadd.txt
