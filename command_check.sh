#!/bin/bash
#BSUB -q priority				#submit to 'priority' queue
#BSUB -n 1						#each  job run on 1 core
#BSUB -W 72:00					#job run 12 hour
#BSUB -o out_check.txt 			#lsf output file
#BSUB -e err_check.txt 			#lsf error file
#BSUB -R 'rusage[mem=8000]'		#use 4GB memory
Rscript './example/glrm_mixing_eval/glrm_mixing_cheker_geweke.R' 
