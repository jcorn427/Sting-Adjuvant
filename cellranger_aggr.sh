#!/bin/sh
#Author John Cornelius 
#Description Cellranger Aggr STING 02 samples 

module load cellranger/8.0.1

cellranger aggr --id=sting_02_full_h5 --csv=/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/cellranger_aggr.csv &
wait