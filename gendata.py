#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Generate data for plots

import os

# Dipole models
models=["ipsat", "ipnonsat", "iim"]
# Gluon distributions
gdists=["dglap"]
#Q^2 values
Q2vals=[0,3,20]

A=197 # Gold
N=200
mint=0
maxt=1.0
num_of_threads=6
bjorkx=0.0001

cmd = "OMP_NUM_THREADS="+str(num_of_threads) + " ./dipole -A " + str(A) \
    + " -N " + str(N) + " -mint " + str(mint)  + " -maxt " + str(maxt) \
    + " -x " + str(bjorkx)
for mode in models:
    for gd in gdists:
        for q in Q2vals:
            filename = mode+"_" + gd + "_" + str(q) + ".txt"
            fullcmd = cmd + " -Q2 " + str(q) + " -dipole " + mode \
                        + " -gdist " + gd + " > data/" + filename
            print (fullcmd)
            os.system(fullcmd)
            os.system("rm -fr tmp; mkdir tmp")
            os.system("sort -n " + filename + " > " + "tmp/" + filename)

print("Done, sorted data can be found under tmp directory")
