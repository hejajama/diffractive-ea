#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Generate data for  d\sigma^A/dt / A*d\sigma^p/dt

import os

# Dipole models
models=["ipsat", "iim", "ipsat_nonsatp"]
# Gluon distributions
gdists=["dglap"]
qvals=[0,10]

wavef="boosted-gaussian"
A=197 # Gold
N=200
t=0.5
num_of_threads=4
gd="dglap"
minx=1e-6
maxx=0.01

for q in qvals:
    cmd = "OMP_NUM_THREADS="+str(num_of_threads) + " ./dipole -A/p_x -A " + str(A) \
        + " -N " + str(N) + " -t " + str(t) \
        + " -minx " + str(minx) + " -maxx " + str(maxx) + " -Q2 " + str(q) \
	+ " -wavef " + wavef 
    for mode in models:
        filename = "data/ap/" + mode + "_x_q" + str(q) + ".txt_"
        fullcmd = cmd + " -dipole " + mode \
                     + " -gdist " + gd + " > " + filename + " "
        print (fullcmd)
        os.system(fullcmd)
        os.system("sort -n " + filename + " > data/ap/" + mode + "_x_q" + str(q) + ".txt")

print("Done")
