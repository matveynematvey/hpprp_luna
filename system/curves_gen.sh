#!/bin/sh

zmin=3.2
zmax=5.1

tmp_script=`tempfile`

echo 'set term png
set yrange['$zmin':'$zmax']
set output '\"$1/curves.png\"'
set grid ytics lw 1 lt 0
set grid xtics lw 1 lt 0
unset key' > $tmp_script

first=true

for x in `ls $1/*.xls`; do

$first && echo "plot '$x' with lines lw 4\\" >> $tmp_script
$first || echo ", '$x' with lines lw 4\\" >> $tmp_script
first=false

done

gnuplot $tmp_script 2>/dev/null

rm $tmp_script

