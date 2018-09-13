#!/bin/bash

ngens=25
nmuts=17

dir=log_files
mkdir $dir

for i in $(seq 1 $ngens); do
	for j in $(seq 1 $nmuts); do
		cp ${i}/${j}/log_antifreeze_score_v7.32_ev.py_1ucs_${i}.${j}.txt ${dir}/
	done
	echo "Done with Gen. "$i
done
