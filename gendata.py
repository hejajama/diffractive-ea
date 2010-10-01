#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Generate data for plots

import os

# Dipole models
models=["ipsat", "ipsat_nonsatp", "ipnonsat", "iim"]
coherentmodels=["ipsat", "ipsat_nonsatp", "ipnonsat", "iim"]
# Gluon distributions
gdists=["dglap"]
#Q^2 values
Q2vals=[0,10]

A=197 # Gold
N=200
mint_coherent=0
mint=0.1
maxt=0.5
num_of_threads=4
bjorkx=0.0001
gd="dglap"

cmd = "OMP_NUM_THREADS="+str(num_of_threads) + " ./dipole -A " + str(A) \
    + " -N " + str(N)  + " -maxt " + str(maxt) \
    + " -x " + str(bjorkx) # + " -scalex"
for mode in models:
    for q in Q2vals:
        filename = "data/q" + str(q) + "/" + mode + ".txt"
        fullcmd = cmd + " -Q2 " + str(q) + " -dipole " + mode \
                     + " -gdist " + gd + " -mint " + str(mint) \
                     + " > " + filename + "_tmp"
        print (fullcmd)
        os.system(fullcmd)
        os.system("sort -n " + filename +"_tmp" + " > data/q" + str(q) + "/" + mode + ".txt")
for mode in coherentmodels:
    for q in Q2vals:
        filename = "data/q" + str(q) + "/" + mode + "_coherent.txt"
        fullcmd = cmd + " -coherent_dt -Q2 " + str(q) + " -dipole " + mode \
                    + " -gdist " + gd + " -mint " + str(mint_coherent) \
                    +  " > " + filename + "_tmp"
        print (fullcmd)
        os.system(fullcmd)
        os.system("sort -n " + filename +"_tmp" + " > data/q" + str(q) + "/" + mode + "_coherent.txt")


print("Done")
