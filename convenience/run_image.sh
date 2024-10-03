#!/bin/bash

arr=(10000 100000)

for i in "${arr[@]}"; do
        for x in {0..10..1}; do
                python make_image_single_wavelength_${x}_${i}.py
        done
done
