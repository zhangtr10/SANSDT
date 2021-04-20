import numpy as np
from subprocess import check_output
import SANS_analysis_t as sans
import cma
import time
import sys
import subprocess
import func
from create_jobs import create_jobs


###variable width of each layer, except the last one and L AL thickness
dim = 18#6 number of layers*3, 3 things need to control: Height, width including Al thickness and Al thickness 

x = [34.92781757, 60.02555013, 40.58570506, 97.03895383, 121.76430538,
     86.56283885, 296.01289259, 187.36359971, 174.88845373, 184.55443451,
     183.28611147, 157.80876657,  54.81904145,  64.46460875,  72.6503356,
     73.9726443,   72.84288907,  63.07784878];#first 6 are heights, next 6 are width and last 6 are 

tolfun = 1e-5##no need to change
popsize=10##number of simlations could be run at the same time
bufsize = 0
iter_cn = 0
h_min = 0.00015##fitting results threshold, as long fitting results are smaller than this value, everything stops.

ferr = open('logerr.log', 'w')
ftemps = open('tmpsol.log', 'w')
logf = open('cma.log', 'w')

#bestever = cma.optimization_tools.BestSolution()
es = cma.CMAEvolutionStrategy( x, 2,{'popsize':popsize,'tolfun':tolfun})## 2 here is the standard deviation of the parameters in x 

while not es.stop():

	cur_fit = np.zeros((1,popsize))##don't change
	solutions =es.ask()##don't change

	func.ftn(dim, solutions)##don't change
	tmp_cnt = 0

	while (tmp_cnt<popsize):
		tmp_cnt = 0
		time.sleep(6000)####waiting for all the jobs finish in one round and then start send round
		for i in range(0,popsize):
			a =func.check_ftn(solutions[i], 'empty_0d_10d_13m.ABS', 'empty_0d_7d_4m.ABS')##'empty_0d_10d_13m.ABS', 'empty_0d_7d_4m.ABS'are experiments files
			
			cur_fit[0][i] = a
			tmp_cnt +=1
		
		logf.write('tmp_cnt '+str(tmp_cnt)+'\n')

	logf.write('succ in cur'+'\n')


	cur_fit = np.array(cur_fit).tolist()[0]

	print(cur_fit, flush=True)##print all the fittings restuls
	
	es.tell(solutions, cur_fit)
    
	tmp_rst = es.result[1]
	ferr.write(str(tmp_rst)+'\n')
	
	tmp_sol = es.result[0]
	print(tmp_sol, flush=True)##print x with the best fittings restuls
	for i in range(0,dim):
	    ftemps.write(str(tmp_sol[i])+' ')
	ftemps.write('\n')

	iter_cn +=1
	logf.write('iter'+str(iter_cn)+'\n')
	
	if(tmp_rst <=h_min ):
	   	break
print (finished)
   
