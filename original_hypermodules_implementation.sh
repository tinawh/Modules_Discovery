#!/usr/bin/env bash

cd /.mounts/labs/reimandlab/private/users/thuang/bin_3/original_hypermodules_implementations/

for file in *.sh; do 
	qsub -cwd -V -l h_vmem=80G $file
done
