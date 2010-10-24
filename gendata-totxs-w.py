#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Generate data for totxs plots

import os

# Dipole models
models=["ipsat", "ipsat-nofactor"]

Q2=0
wavef=["boosted-gaussian", "gaus-lc"]
#wavef="gaus-lc"
A=1
N=200
mint_coherent=0
minw=20
maxw=300
num_of_threads=1
num_of_processes=4
gd="dglap"

process=0

cmd = "OMP_NUM_THREADS="+str(num_of_threads) + " ./dipole -A " + str(A) \
    + " -N " + str(N) + " -totxs_w -minW " + str(minw) \
    + " -maxW " + str(maxw) + " -Q2 " + str(Q2) + " " 
for mode in models:
    for wave in wavef:
        filename = "data/totxs/" + mode + "_" + wave + "_w.txt"
        fullcmd = cmd + "-dipole " + mode \
                     + " -wavef " + wave + " > " + filename + "_tmp"
        process=process+1
	if (process<num_of_processes):		# Run in background
		fullcmd = fullcmd + " &"
	else:
		process=0
	print (fullcmd)
        #os.system(fullcmd)
        

for mode in models:
	for wave in wavef:
		filename = "data/totxs/" + mode + "_" + wave + "_w.txt"
		sortcmd = "sort -n " + filename +"_tmp" + " > data/totxs/" + mode + "_" + wave + "_w.txt"
		print (sortcmd)
		os.system(sortcmd)

print("Done")
