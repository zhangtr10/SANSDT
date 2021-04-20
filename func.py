import cma
import numpy as np 
import sys
import os.path
import subprocess
import SANS_analysis_t as sans
from create_jobs import create_jobs
import pandas as pd
import time

def ftn_results(file1,file2,file3,file4):
##load experiment data file1_13m and file3_4m of different deterctor distance
  empty_exp_13m = sans.experiment()
  empty_exp_13m.read_1D_sector(file1)

  empty_exp_4m = sans.experiment()
  empty_exp_4m.read_1D_sector(file3)
##load data from models
  model_data = sans.model(file2)
  model_data.sector_average_full(10)

  model_data_4m = sans.model(file4, detector_distance=4000)
  model_data_4m.sector_average_full(10)
##extract intensity as a function of q from experiments 
  x_exp, y_exp = empty_exp_13m.out_put_data_exp(trim_lead=8,shift_up=(1./3.5E-4))
  x_exp=np.log10(x_exp)
  y_exp=np.log10(y_exp)

  x_exp_4m, y_exp_4m = empty_exp_4m.out_put_data_exp(trim_lead=4,shift_up=(1./3.5E-4))
  x_exp_4m=np.log10(x_exp_4m)
  y_exp_4m=np.log10(y_exp_4m)
##extract intensity as a function of q from models 
  x_sim, y_sim =model_data.out_put_data(trim_lead=10, trim_tail=49, shift_up=(1./3.5E-4)*8.5E3)

  x_sim=np.log10(x_sim)
  y_sim=np.log10(y_sim)

  x_sim_4m, y_sim_4m =model_data_4m.out_put_data(trim_lead=7, trim_tail=59, shift_up=(1./3.5E-4)*7.5E2)

  x_sim_4m =np.log10(x_sim_4m)
  y_sim_4m =np.log10(y_sim_4m)

  y_interp_exp = np.interp(x_sim, x_exp,y_exp )
  y_interp_exp_4m = np.interp(x_sim_4m, x_exp_4m,y_exp_4m )
##fitting resutls, criteria to confirm the goodness of the fittings. The smaller the better(<0.001)
  fit_1=np.trapz((y_sim-y_interp_exp)*(y_sim-y_interp_exp), x=x_sim)/np.trapz(y_interp_exp*y_interp_exp, x=x_sim)
  fit_2=np.trapz((y_sim_4m-y_interp_exp_4m)*(y_sim_4m-y_interp_exp_4m), x=x_sim_4m)/np.trapz(y_interp_exp_4m*y_interp_exp_4m, x=x_sim_4m)
  print(fit_1,fit_2,flush=True)##print fitting resutls for different models
  return (fit_1+fit_2)/2

def get_directory_name(solution_i):

  return '_'.join('{:.2f}'.format(solution_ii) for solution_ii in solution_i)

def get_directory_name_4m(solution_i):

  return '__'.join('{:.2f}'.format(solution_ii) for solution_ii in solution_i)

def ftn(dim,x):
  #import pdb; pdb.set_trace()
  ##initilize the jobs to run on the cluster. run 13m and 4m jobs together since we want to fit both at the same time
  input_files = ['input_al', 'input_al_4m','submit_PFT','submit_PFT_4m']
  
  job_names = [get_directory_name(x_i) for x_i in x]
  job_names_4m = [get_directory_name_4m(x_i) for x_i in x]

  heights = [' '.join(list(map(str, abs(x_i[:6])))) for x_i in x]
  widths = [' '.join(list(map(str, abs(x_i[6:12])))) for x_i in x]
  #al_thicknesses = [abs(x_i[8]) for x_i in x]
  al_thicknesses = [' '.join(list(map(str, abs(x_i[12:])))) for x_i in x]
  df = pd.DataFrame({'JOB_NAME': job_names, 'WIDTHS': widths, 'HEIGHTS': heights, 'AL_THICKNESS': al_thicknesses})
  create_jobs(input_files, df, sub_file='submit_PFT')
  

  df_4m = pd.DataFrame({'JOB_NAME': job_names_4m, 'WIDTHS': widths, 'HEIGHTS': heights, 'AL_THICKNESS': al_thicknesses})
  create_jobs(input_files, df_4m, sub_file='submit_PFT_4m')#create folders and submit jobs meantime.

  # file=str(int(x[0]))
  # for i in range(1,dim):
  #   file+='_'+str(int(x[i]))
  #   file+='_'+'test'

  # out = subprocess.check_output(["pwd"])
  # out2 = out[0]
  # file2 = write_input(dim, x)

  # for  i in range(1,len(out)-1):
  #   out2 += out[i]
   
  # file=out2+'/'+file
  # sub_t = ["./sub_pmf.sh",file]
   
  # file3=file2+'_output'
  # sub_t.append(file2)
  # sub_t.append(file3)
   #print sub_t
   
  #subprocess.call(sub_t)

def write_input(dim, x):
  file=str(int(x[0]))
  for i in range(1,dim):
    file+='_'+str(int(x[i]))
  file+='_'+'input'
  text=open(file,"w")
  text.write("%d\n"% (9))
  text.write("%d\n"% (10000))
  text.write("%d\n"% (25))
  text.write("%d\n"% (360))
  text.write("%d\n"% (400))
  text.write("%d\n"% (190))
  text.write("%f\n"% (abs(x[int(dim/2-1)])))
  text.write("%d\n"% (13000))

  for i in range(int(dim/2)):
    text.write("%f "% (abs(x[i])))
  text.write("\n")
  for i in range(int(dim/2)):
    text.write("%f "% (abs(x[i+int(dim/2)])))  	
  return file
 


def check_ftn(x, file1, file2):
  sim_dir = get_directory_name(x)
  sim_dir_4m = get_directory_name_4m(x)

  ## check the job status: finish or not
  while not (os.path.isfile(os.path.join(sim_dir, 'output'))):
    time.sleep(2000)

  while not (os.path.isfile(os.path.join(sim_dir_4m, 'output_4m'))):
    time.sleep(2000) 
## return fitting rueslts
  return ftn_results(file1, os.path.join(sim_dir, 'output'),file2,os.path.join(sim_dir_4m, 'output_4m'))
  

def remove(dim,x):
    file='./'+str(int(x[0]))
    for i in range(1,dim):
        file+='_'+str(int(x[i]))
    file+='_'+'test'
    subprocess.call(["rm","-r",file])
















