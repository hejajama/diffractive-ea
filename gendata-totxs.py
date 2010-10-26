#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Generate data for totxs plots

import os

# Dipole models
models=["ipsat", "ipnonsat", "iim", "ipsat-nofactor"]

# Gluon distributions
gdists=["dglap"]

W=90
wavef="boosted-gaussian"
#wavef="gaus-lc"
A=1
N=200
mint_coherent=0
minq2=0
maxq2=90
num_of_threads=4
bjorkx=0.0001
gd="dglap"
processes=3
running=0

cmd = "OMP_NUM_THREADS="+str(num_of_threads) + " ./dipole -A " + str(A) \
    + " -N " + str(N) + " -totxs_q2 -minQ2 " + str(minq2) \
    + " -maxQ2 " + str(maxq2) + " -W " + str(W) + " " \
    + "-wavef " + wavef + " " 
for mode in models:
    filename = "data/totxs/" + mode + ".txt"
    fullcmd = cmd + "-dipole " + mode \
                     + " -gdist " + gd + " > " + filename + "_tmp"
    if (running<processes):
         fullcmd = fullcmd + " &"
         running=running+1
    else:
         running=0
    print (fullcmd)
    os.system(fullcmd)

for mode in models:
    filename = "data/totxs/" + mode + ".txt"
    os.system("sort -n " + filename +"_tmp" + " > data/totxs/" + mode + ".txt")

print("Done")
