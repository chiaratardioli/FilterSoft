#!/bin/bash

usage () {
        echo " Usage: chain_order.sh <no. of alias for each object - 1>"
#        echo " comments ..."
        exit 1
}
if [ $# -lt 1 ]; then
        echo 'no number supplied'
        usage
fi
if [ $# -gt 865 ]; then
        echo 'number too large'
        usage
fi

cd datafile;
make clean;

cd ../crea_cat;
make clean
echo $1|./lenz_cat.x;

if [ $# -gt 0 ]; then
    let "norbs = $((($1*2+1)*864))";
else 
    let "norbs = 864";
fi
echo 'no. of orbit to process:' $norbs;

./genera_inputs.pl $norbs; #use file allorbits.fla given by genera.x

cd ../order;
mv input.in input.in.aux;

for ((i=1;i<=$norbs;i++)); do
mv ../crea_cat/input.in.$i .;
mv input.in.$i input.in;
./order>output;
./insertline.py;
mv tmp_results ../datafile/results_kepl_$i.dat;
done;
mv input.in.aux input.in;

cd .. ;
#echo $norbs | ./evolution_file.x; (old)
echo $norbs | ./ephem_table.x;
