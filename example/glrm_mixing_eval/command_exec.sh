#!/bin/bash
#BSUB -q short					#submit to "short" queue
#BSUB -n 1						#each  job run on 1 core
#BSUB -W 12:00					#job run 12 hour
#BSUB -J jobArray[1-900]		#job array list goes 1,2,3...n
#BSUB -o out_%I.txt 			#lsf output file
#BSUB -e err_%I.txt 			#lsf error file
#BSUB -R 'rusage[mem=4096]'		#use 4GB memory
Rscript './example/glrm_mixing_cook/glrm_mixing_master.R' $LSB_JOBINDEX