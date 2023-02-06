#!/bin/bash


rm -f .jobs.txt
rm -r ~/Analysis/logfiles
rm -r ~/Analysis/outputs

mkdir ~/Analysis/logfiles
mkdir ~/Analysis/outputs

cp ~/ranscripts/.jobs.txt ~/Analysis/.jobs.txt

xargs.sh .jobs.txt
