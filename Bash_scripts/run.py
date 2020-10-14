# simulations on residuals

import pyreadr
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import argparse
import time
from datetime import datetime
import multiprocessing as mp
import os


def get_genotype_vector(freq, N=369153):
    genotypes = np.random.choice([0,1,2], N, p=[freq**2, 2*freq*(1-freq), (1-freq)**2])
    return genotypes


def run_n_sim(freq, N, data):
    add_list = []
    log_list = []
    for i in range(N):
        data['gen'] = get_genotype_vector(freq)
        add_p = smf.ols(formula='add ~ gen', data=data).fit().pvalues['gen']
        log_p = smf.ols(formula='log ~ gen', data=data).fit().pvalues['gen']
        add_list.append(add_p)
        log_list.append(log_p)	
    return pd.DataFrame({'add': add_list, 'log': log_list})


def get_simulation_results(freq, number_of_sim, data, n_threads=2):
    result_list = []

    def log_result(result):
        # This is called whenever foo_pool(i) returns a result.
        # result_list is modified only by the main process, not the pool workers.
        result_list.append(result)

    def apply_async_with_callback():
        pool = mp.Pool(n_threads)
        for i in range(n_threads):
            pool.apply_async(run_n_sim,
                             args=(freq, number_of_sim//n_threads, data +1),
                             callback=log_result)
        pool.close()
        pool.join()

    start_time = time.time()
    apply_async_with_callback()
    print(f"{time.time() - start_time} seconds to do {number_of_sim} simulations, frequency {freq}")
    output = pd.concat(result_list)
    return output

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-nsim', dest='nsim', type=int, help="number of simulations for each frequency", required=True)
	parser.add_argument('-nthreads', dest='n_threads', type=int, help="number of threasds", required=True)
	parser.add_argument('-freq_list', nargs='+', help='frequencies list', required=True)
	parser.add_argument('-outpath', dest='outpath', type=str, default='./sim_res',
		help="path to store output", required=False)
	parser.add_argument('-rdatapath', dest='rdatapath', type=str, default='./residuals.RData',
		help="number of simulations for each frequency", required=False)


	args = parser.parse_args()

	data = pyreadr.read_r(args.rdatapath)['resid']

	if not os.path.exists(args.outpath):
		os.makedirs(args.outpath)	

	for freq in args.freq_list:
		df = get_simulation_results(freq=float(freq), number_of_sim=args.nsim, data=data, n_threads=args.n_threads)
		date_str = datetime.now().strftime("%m_%d_%H_%M")
		tsv_name = f'res_{freq}_{date_str}.tsv'
		df.to_csv(os.path.join(args.outpath, tsv_name), sep='\t', index=False)
