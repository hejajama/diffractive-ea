#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Generate data for  d\sigma^A/dt / A*d\sigma^p/dt

import os

# Dipole models
models=["ipsat", "ipsat_nonsatp", "iim"]
# Gluon distributions
gdists=["dglap"]
xvals=[0.0001, 0.01]

A=197 # Gold
N=100
t=0.5
num_of_threads=1
#bjorkx=0.0001
gd="dglap"
wavef="boosted-gaussian"
maxq2=100
minq2=0.01

for bjorkx in xvals:
    cmd = "OMP_NUM_THREADS="+str(num_of_threads) + " ./dipole -A/p -A " + str(A) \
        + " -N " + str(N) + " -t " + str(t) \
        + " -x " + str(bjorkx) + " -maxQ2 " + str(maxq2) + " -minQ2 " + str(minq2) \
	+ " -wavef " + wavef
    for mode in models:
        filename = "data/ap/" + mode + "_x" + str(bjorkx) + ".txt_"
        fullcmd = cmd + " -dipole " + mode \
                     + " -gdist " + gd + " > " + filename 
        print (fullcmd)
        os.system(fullcmd)
        os.system("sort -n " + filename + " > data/ap/" + mode + "_x" + str(bjorkx) + ".txt")

print("Done")
