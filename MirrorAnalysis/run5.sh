#!/bin/bash


if [ -z ${1+x} ]; then echo "give a hld file name as argument"; exit; fi
./MirrorAnalysis.x -t hld -f $1 -p conf_trb3.xml -c TOTConfigRun5After.root -u userParams5Run.json -l detectorSetupRun2_3_4_5.json -i 2 -o ~/JPET_Data/ -r 1 400000
