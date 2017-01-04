cd Project/2.BayesOpt/

module load stats/R/3.3.1
module load dev/atlas/3.10.2 
module load dev/blas

# remove archived files
find . -name "Archive" -type d -ls -exec rm -rv {} +


# examine disk space, don't exceed 100G
du -h ./

# delete result
find ./result/ -name "*.RData" -print0 | xargs -0 rm

# clean up and execute
#rm err*.* out*.*
#rm ./result/mixing_stat/*.*
bkill -b 0

bsub < ./command_exec.sh 

bsub -q priority -W 72:00 -o out_sum.txt -e err_sum.txt 
	 -R 'rusage[mem=4096]' -N Rscript "./example/glrm_mixing_cook/glrm_mixing_cheker.R"; 
