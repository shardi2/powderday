#!/bin/bash

#path  = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameters_master_"

#model_files = "/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameters_model_*"

for i in {0..10..1}
do
	python sed_compare_photons_10mil_halo_${i}.py

done	
