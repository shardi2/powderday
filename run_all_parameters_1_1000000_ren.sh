#!/bin/bash

for i in {0..29..1}
do
	python pd_front_end.py /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_renaissance/parameter_files_1 parameters_master_${i}_1000000 parameters_model_${i}_1000000 
done	
