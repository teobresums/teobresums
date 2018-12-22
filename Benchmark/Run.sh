#!/bin/bash
# Script to run all the parfiles

for p in `ls ./par/*.par`; do
    name=$(basename "$p")
    name="${name%.*}"
    echo "$name"
    ./TEOBResumS.x $p &> $name.out
    sed -i "s/output_dir = .*/output_dir = \"$name\"/" $p
    mv $name/ data/
done

