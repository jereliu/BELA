#!/bin/bash
#BSUB -q short			#submit to 'priority' queue
#BSUB -n 1					#each  job run on 1 core
#BSUB -W 12:00				#job run 12 hour
#BSUB -J jobArray[1-480]		#job array list goes 1,2,3...n
#BSUB -o './log/out_check_%I.txt' 		#lsf output file
#BSUB -e './log/err_check_%I.txt' 		#lsf error file
#BSUB -R 'rusage[mem=8000]'	#use 4GB memory
Rscript './example/glrm_mixing_eval/glrm_mixing_cheker_geweke.R' $LSB_JOBINDEX