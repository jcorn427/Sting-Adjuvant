#!/bin/sh
#Author John Cornelius 
#Description Cellranger Multi STING 02 samples 

module load cellranger/7.1.0

cellranger multi  --id=sting_02_part_2_set2_C --csv=/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/multi_config_sting-adjuvant_02_part_2_set2_C.csv &
wait