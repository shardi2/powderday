#!/bin/bash

#path  = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameters_master_"

#model_files = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameters_model_*"

for i in {725..754..1}
do
	python pd_front_end.py /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_renaissance/parameter_files_26 parameters_master_${i}_1000000 parameters_master_neb_${i}_1000000 parameters_model_${i}_1000000 
done	
