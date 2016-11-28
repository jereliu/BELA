Functions to evaluate on convergence using method from Cook Gelman and Rubin (2006)

1. glrm_mixing_master: 
	* select random seed for each worker

2. glrm_mixing_worker: 
	* take seed, generate data, run model
	* calculate test statistics for each iteration
	* return iteration-indexed list of test statistics

3. glrm_mixing_checker:
	* summarize result into array of [iter, stat]
	* calculate sup|F'(theta) - F(theta)| for each iteration.