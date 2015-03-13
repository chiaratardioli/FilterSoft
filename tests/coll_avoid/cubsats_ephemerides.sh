#!/bin/bash

#usage () {
#        echo " Usage: chain_order.sh <norb>"  # echo " comments ..."
#        exit 1
#}
#if [ $# -lt 1 ]; then
#        echo 'no number supplied'
#        usage
#fi
#if [ $# -gt 865 ]; then
#        echo 'number too large'
#        usage
#fi

cd datafile;
make clean;

cd ../cubesat;
rm -f allorbitsSD.*;
#rm -f file.txt;
#ln -s cubesat.txt file.txt;

# To be set before
#cat cubesat.txt molniya.txt >> file.fla

./TLE2KEP.x
mv file_oe.fla allorbitsSD.fla;
NUMLIN=$(cat allorbitsSD.fla | wc -l ); # no. of lines
echo "No. of objects: $NUMLIN";
./genera_inputs.pl $NUMLIN; #use file allorbitsSD.fla given by the space debris

cd ../order;
mv input.in input.in.aux;

for ((i=1;i<=$NUMLIN;i++)); do
mv ../cubesat/input.in.$i .;
mv input.in.$i input.in;
./order > output;
./insertline.py;
mv tmp_results ../datafile/results_kepl_$i.dat;
#mv results_kepl.dat ../datafile/results_kepl_$i.dat;
#mv tmp_results ../../coll_avoid/datafile/results_kepl_$i.dat;
done;
mv input.in.aux input.in;

#cd ../orb41/tests/coll_avoid;
cd ..;
echo $NUMLIN | ./evolution_file.x
