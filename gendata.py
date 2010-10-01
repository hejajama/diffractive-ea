#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Generate data for plots

import os

# Dipole models
models=["ipsat", "ipnonsat", "iim", "ipsat -coherent_dt", "ipnonsat -coherent_dt"]
# Gluon distributions
gdists=["dglap"]
#Q^2 values
Q2vals=[0,10]

A=197 # Gold
N=200
mint=0
maxt=0.5
num_of_threads=4
bjorkx=0.0001
gd="dglap"

cmd = "OMP_NUM_THREADS="+str(num_of_threads) + " ./dipole -A " + str(A) \
    + " -N " + str(N) + " -mint " + str(mint)  + " -maxt " + str(maxt) \
    + " -x " + str(bjorkx) # + " -scalex"
for mode in models:
    for q in Q2vals:
        filename = "data/q" + str(q) + "/" + mode + ".txt"
        fullcmd = cmd + " -Q2 " + str(q) + " -dipole " + mode \
                     + " -gdist " + gd + " > " + filename + "_tmp"
        print (fullcmd)
        os.system(fullcmd)
        os.system("sort -n " + filename +"_tmp" + " > data/q" + str(q) + "/" + mode + ".txt")
print("Done")
