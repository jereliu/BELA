source("./example/glrm_mixing_eval/glrm_mixing_worker.R")

#### 0. Read-In Call Number ####

# read config
print("Mr Handy: How may be of service, master?")
print("Mr Handy: Oh a 'argument', how wonderful..")
args <- commandArgs(trailingOnly = TRUE)
job_idx <- as.numeric(args)
print(sprintf("Mr Handy: You see the job index '%d'..", job_idx))
print("Mr Handy: and who gets to read all this mumble jumble? Me, that's who...")

#### 1. Config Generation/Read In ====
n_data <- 1
n_rep <- 1e3
n_run_per_worker <- 5

cfig_file <- "cfigList.csv"
if (!file.exists(cfig_file)){
  set.seed(100)
  
  # create data list
  cfig_list <-
    expand.grid(
      data_seed = sample(n_data*1e3, n_data),
      trial_seed = sample(n_rep*1e3, n_rep),
      K_true = c(15, 25), 
      K_model = c(5, 10, 15, 20, 25, 30),
      SNR = 100, LAMBDA = 10,
      FAMILY = c("gaussian", "poisson"),
      SAMPLR = c("hmc_stan", "vi_stan")
    )
  
  write.csv(cfig_list, file = cfig_file, 
            row.names = FALSE)
  
  # produce execution file
  n_exec <- 
    ceiling(nrow(cfig_list) / (n_run_per_worker * 1e3))
  for (exec_id in 1:n_exec) {
    n_rask_rest <- 
      min(nrow(cfig_list)/n_run_per_worker - 
            (exec_id-1)*1e3, 1e3)
    string <- 
      paste0(
        "#!/bin/bash
  #BSUB -q short					#submit to 'short' queue
  #BSUB -n 1						#each  job run on 1 core
  #BSUB -W 12:00					#job run 12 hour
  #BSUB -J jobArray[", (exec_id-1) * 1e3 + 1, 
        "-", (exec_id-1) * 1e3 + n_rask_rest, 
        "]		#job array list goes 1,2,3...n
  #BSUB -o './log/out_%I.txt' 			#lsf output file
  #BSUB -e './log/err_%I.txt' 			#lsf error file
  #BSUB -R 'rusage[mem=4096]'		#use 4GB memory
  Rscript './example/glrm_mixing_eval/glrm_mixing_master.R' $LSB_JOBINDEX"
      )
    write(string, file = paste0("command_exec", exec_id, ".sh"))
  }
}
cfig_list <- read.csv(cfig_file)

#### 2. Execution ====
config_idx_list <-
  ((job_idx - 1) * n_run_per_worker + 1):
  (job_idx * n_run_per_worker)

for (config_idx in config_idx_list){
  cfig <- cfig_list[config_idx, ]
  
  #### 1. run worker ####
  N <- 50
  P <- 100
  
  record_freq <- 
    ifelse(cfig$SAMPLR == "vi_stan", 5, 1)
  iter_max <- 
    ifelse(cfig$SAMPLR == "vi_stan", 1e3, 5e2)
  
  status <- 
    glrm_worker(
      data_seed = cfig$data_seed, 
      trial_seed = cfig$trial_seed,
      n = N, p = P, 
      K_true = cfig$K_true,
      K_model = cfig$K_model,
      SNR = cfig$SNR,
      LAMBDA = cfig$LAMBDA,
      FAMILY = cfig$FAMILY,
      SAMPLR = cfig$SAMPLR,
      iter_max = rep(iter_max, 2),
      record_freq = record_freq,
      parm_updt = c("U", "V")
    )
}

print("Mr Handy: That it? Ah, fine then.")