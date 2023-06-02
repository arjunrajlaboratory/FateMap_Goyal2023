#! /bin/bash

module purge
module load singularity
module load nextflow

nextflow run nf-core/sarek -r 3.0 -work-dir ./sarek/work -params-file nf-params.json 