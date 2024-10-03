#!/bin/bash

arr=(10000 100000)

for i in "${arr[@]}"; do
	for x in {0..10..1}; do
		python sed_plot_${x}_${i}.py
	done
done
