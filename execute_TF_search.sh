#!/bin/sh
##$BATCH -l h_vmem=200G
#$BATCH --mem=80G
#purge previous loads
module purge
#start loading the path
module load modulepath
#load the modules
module load Python/3.8.6-GCCcore-10.2.0

python get_TF.py

