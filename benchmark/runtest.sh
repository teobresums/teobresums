#!/bin/bash
# Simple script to benchmark the code

SRCPATH=/home/bernuzzi/Codes/teobresums/C++/
MAKEFILE=/home/bernuzzi/Codes/teobresums/C++/Makefile.TEOBResumS

#
#
#
NOW=$(date +"%Y%m%d%H%M%S")
TOP=$(pwd)
GIT=$(git rev-list HEAD | wc -l)
NAME="bench_${GIT}_${NOW}"

#
#
#
mkdir -pv $NAME
cd $NAME
WORKDIR=$(pwd)
echo $(pwd)
cp -r $SRCPATH/src .
cp $MAKEFILE Makefile
make 
tar zcvf src.tgz Makefile src/
rm -rvf src/ obj/ Makefile
git log --format="%H" -n 1 > git_status.txt

#
#
#
for p in `ls ../parfiles/*.par`; do
    test=$(basename "$p")
    test="${test%.*}"
    echo "$test"
    ./TEOBResumS -p $p > $test.out
    mkdir -p data/$test/
    mv waveform.dat data/$test/
done
