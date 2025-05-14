#!/bin/bash

for i in {1421..1450..1}
do
	python pd_front_end.py /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_renaissance/parameter_files_50 parameters_master_${i}_1000000 parameters_model_${i}_1000000 
done	
