#!/bin/bash


rm -f .jobs1.txt
rm -r ~/Analysis/logfiles_1ak8
rm -r ~/Analysis/output_1ak8

mkdir ~/Analysis/logfiles_1ak8
mkdir ~/Analysis/output_1k8

cp ~/Analysis/tools/.jobs1.txt ~/Analysis/.jobs1.txt

xargs.sh .jobs1.txt
