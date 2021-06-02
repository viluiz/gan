from __future__ import print_function 
import numpy as np
import sys
import os
import random
import time


''' This file with create a number of directories. It will create input files
    with values of R0 inside. It will then run the fortran virus modelling code
    and then compress the results.'''

#########################################################################
# variables that must be set
nruns = 200 # number of simulations
nrandom = 2 # 2 or 10 # number of variables set randomnly per simulation
random_sample_method = 'normal' # 'normal' or 'uniform'
mean = 10; std = 4; cut_off = 0.2 # only used for normal

#########################################################################

if nrandom !=2 and nrandom !=10:
    print ('nrandom can only be 2 or 10! Aborting...') 
    sys.exit()

# two hardwired cases
if nrandom ==2:
    random_names = ['R0_home','R0_mobile', ] # names of the random variables
    executable = '../diffusion3d_virus > /dev/null & ' 

elif nrandom == 10: 
    random_names = ['R0_home','R0_office', 'R0_school', 'R0_hospital', 'R0_pedest', 'R0_diff1', 'R0_diff2',  'R0_shops', 'R0_vehc1to5', 'R0_park'] # 
    executable = '../diffusion3d_virus_town_largg > /dev/null & '

#########################################################################
t_compute = 0.0
for irun in range(nruns):

    jrun = irun + 1

    print('simulation number', jrun)

    # move to subdirectory (create if necessary)
    directory = 'results_run_' + str(jrun) 
    #print('pwd', os.getcwd())
    if not os.path.exists(directory):
        os.mkdir(directory)
    os.chdir(directory)
    #
    #print('pwd', os.getcwd())


    # set the R0s by selecting random numbers and scaling
    if random_sample_method == 'uniform':
        random_values = []
        for k in range(nrandom): 
            random_values.append(random.random()) 
        random_variables = [0.2+19.8*r  for r in random_values]
        # random_variables is a list

    elif random_sample_method == 'normal':
 
        random_variables = np.random.normal(mean, std, nrandom)
        random_variables[random_variables<cut_off] = cut_off
        # random_variables is a np array

    print('R0_names:')
    print(random_names)
    print('R0_variables')
    print(random_variables)

    # write input file
    # header
    input_file = open('r0-' + str(nrandom) + 'values.csv', 'w')
    line = '# ' + ' '.join(random_names) +'\n'
    input_file.write(line)
    # random variables
    string_list = [str(i) for i in random_variables]
    line = ', '.join(string_list) +'\n'
    input_file.write(line)
    input_file.close()

    # run the code
    t0 = time.time()
    print('command',executable)
    os.system(executable)
    t1 = time.time() 
    t_compute = t_compute + t1 - t0
    print('')

    # compress data
    #os.system('gzip run/group-output-time' + str(irun) + '.csv') #! factor of 12 reduction in disc space.

    os.chdir('../')
    #print('pwd', os.getcwd())

print('Simulations have finished.' , nruns , 'simulations have run in', t_compute , 'seconds.' )
f = open('log.txt',"a")
f.write("Time taken for  %d  simulations to run is %g seconds.\n"  % (nruns, t_compute )  )
f.close()

