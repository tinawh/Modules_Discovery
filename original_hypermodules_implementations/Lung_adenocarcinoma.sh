#!/bin/bash
module load java/1.6.0_21
java -jar /.mounts/labs/reimandlab/private/users/thuang/HyperModules_1.0.2_CMD/HyperModulesCMD-1.0.2.jar -n /.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv -s /.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/mutation_data/mut_Lung_adenocarcinoma_hugo.csv -c /.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/clinical_data/Lung_adenocarcinoma_clin.csv -S 0 -t logrank > /.mounts/labs/reimandlab/private/users/thuang/data_3/original_hypermodules_implementations/Lung_adenocarcinoma.txt
