#!/bin/bash

#path  = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_3000000/parameters_master_"

#model_files = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_3000000/parameters_model_*"

for i in {0..10..1}
do
	python pd_front_end.py /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_3000000  parameters_master_${i}_3000000 parameters_model_${i}_3000000
done	
