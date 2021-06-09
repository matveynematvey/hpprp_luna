#!/bin/sh

zmin=1.2
zmax=9.7

for x in `ls $1/*.xls`; do
	echo 'set term png
set yrange['$zmin':'$zmax']
set output '\"$x.png\"'
set grid ytics lw 1 lt 0
set grid xtics lw 1 lt 0
unset key

plot '\"$x\"' with lines lw 4' | gnuplot - 2>/dev/null
done

mencoder mf://$1/*.png -mf fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $2 >/dev/null 2>/dev/null


tmp_script=`tempfile`

echo 'set term png
set yrange['$zmin':'$zmax']
set output '\"$1/all.png\"'
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

