#!/usr/bin/env bash

cd /.mounts/labs/reimandlab/private/users/thuang/bin_3/hypermodules_implementations/

for f in ./*/*; do qsub -cwd -V -l h_vmem=80G $f; done
