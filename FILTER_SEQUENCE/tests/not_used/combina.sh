#!/bin/bash

usage () {
        echo " Usage: combina.sh <n1> <n2)"
#        echo " comments ..."
        exit 1
}
if [ $# -lt 1 ]; then
        echo 'no numbers supplied'
        usage
fi

rm -f results_kepl_1.dat;
rm -f results_kepl_2.dat;

ln -s datafile/results_kepl_$1.dat results_kepl_1.dat;
ln -s datafile/results_kepl_$2.dat results_kepl_2.dat;
