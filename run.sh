#!/bin/bash


rm -f .jobs.txt
rm -r ~/Analysis/logfiles
rm -r ~/Analysis/output

mkdir ~/Analysis/logfiles
mkdir ~/Analysis/output

cp ~/Analysis/tools/.jobs.txt ~/Analysis/.jobs.txt

xargs.sh .jobs.txt
