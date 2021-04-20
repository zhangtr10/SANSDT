# SANS_DT
Instruction to use optimization:
1. Install c++ Eigen library.
2. For my study, I have two sets of experiments data: one with detector distance of 13m and one with detector distance of 4m.  I have to run fit these two sets of data at the same time. I need to compile two c++ code with command: "g++ -I/~/eigen/ -std=c++11 -O3 neutron_multilayer_input_al.cpp -o a.out" and "g++ -I/~/eigen/ -std=c++11 -O3 neutron_multilayer_input_al_4000.cpp -o b.out"
3. Two c++ codes could read input_al and input_al_am separately and most of the input parameters could be adjusted in these two files.
4. After compiling these two codes, the Covariance Matrix Adaptation Evolution Strategy (CMA-ES) optimization algorithm could be ran in python. You need to install cma package in python  first and a function called create_jobs which can be found on https://github.com/benlindsay/create_jobs . 
For more information about CMA-ES, please refer to https://github.com/CMA-ES/pycma
5. In the python main.py, you could adjust the geometry of the templates for the simulation. And this python program could submit c++ jobs automatically.
If you have any questions regarding the codes, please feel free contact me with email : tianren@seas.upenn.edu
