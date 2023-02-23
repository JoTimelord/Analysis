#!/bin/bash


# rm -f .jobs2.txt
rm -r ~/Analysis/logfiles
rm -r ~/Analysis/output2

mkdir ~/Analysis/logfiles
mkdir ~/Analysis/output2

# cp ~/Analysis/tools/.jobs2.txt ~/Analysis/.jobs2.txt

xargs.sh .jobs3.txt
