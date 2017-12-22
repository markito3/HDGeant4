#!/bin/bash

for i in {-100..100..50}
do
    for j in {-100..100..50}
    do
	echo $i $j
	#	root -l -q -b loadlib.C drawHitPats1.C\($i-5, $i+5, $j-5, $j+5\)
	root -l -q -b loadlib.C drawHitPats1.C\($i, $i+10, $j, $j+10\);
    done
done
	 
echo "All done"
