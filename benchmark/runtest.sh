#!/bin/bash
# Script to run all the parfiles

for p in `ls ./parfiles/*.par`; do
    name=$(basename "$p")
    name="${name%.*}"
    echo "$name"
    ./TEOBResumS.x -p $p > $name.out
    mkdir -p data/$name/
    mv waveform.dat data/$name/
done

