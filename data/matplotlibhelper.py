#!/usr/bin/python
# -*- coding: UTF-8 -*- 

# Reads data from a file in format
# x y

import sys
sys.path.append('/home/hejajama/lib/python');

from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

# 
dashes = ['--', #    : dashed line
          '-', #     : solid line
          '-.', #   : dash-dot line
          ':', #    : dotted line
           '-']


# Possibility to use , as a desimal separator
import locale
locale.setlocale(locale.LC_ALL,"fi_FI")
fnx = lambda x : locale.format("%.1f", x)
fnx2 = lambda x : locale.format("%.2f",x)  
fnx3 = lambda x : locale.format("%.3f",x) 

# Read data and multiply all y values by factor m
def readfile(file, xlist, ylist,m=1.0): 
	f = open(file,"r")
	lines=f.readlines()
	n=len(lines)

	for i in range(n):
		s=lines[i].split()
		if (len(s)==2):
			xlist.append(float(s[0]))
			ylist.append(float(s[1])*m)
	f.close()


