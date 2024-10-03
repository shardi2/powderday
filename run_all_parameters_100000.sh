#!/bin/bash

#path  = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameters_master_"

#model_files = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameters_model_*"

for i in {0..10..1}
do
	python pd_front_end.py /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/ parameters_master_${i}_100000 parameters_model_${i}_100000
done	
