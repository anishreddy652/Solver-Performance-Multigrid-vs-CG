# Solver-Performance-Multigrid-vs-CG



In this project we compare the performances of solvers based on the metrics of execution time and number of iterations taken to convergerve below a certain tolerance limit.
We compare the performances of the following solvers.
1) Multigrid Solver parallelized with OpenACC
2) Red-Black Solver parallelized with OpenACC
3) Conjugate Gradient Solver. 




# Instructions to run parallelized scripts:

First, SSH into the SCC

Get into queue and reserve a GPU in an interactive session
qrsh -l h_rt=01:45:00 -pe omp 1 -P paralg -l gpus=1.0 -l gpu_c=6.0


Load in PGI Compiler & gcc
module load pgi
module load gcc

Compile PGI code multigrid_OpenAcc.cpp
pgc++ -std=c++11 -Minfo=accel -acc -ta=tesla multigrid_OpenAcc.cpp -o runmg

Run
./runmg

Compare with Serial

Compile
g++ -O2 multigrid.cpp -o serial

Run
./serial


Conjugate Gradient can be run from the makefile within the respective directory.


