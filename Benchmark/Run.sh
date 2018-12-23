#!/bin/bash
# Script to run all the parfiles

for p in `ls *.par`; do
    name=$(basename "$p")
    name="${name%.*}"
    echo "$name"
    (time ./TEOBResumS.x $p) &> $name.out
    sed -i "s/output_dir = .*/output_dir = \"$name\"/" $p
done

