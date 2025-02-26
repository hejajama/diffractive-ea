#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Plot dsigma/dt

import sys
sys.path.append("../")
from matplotlibhelper import *


def AddJpsiMass(xlist):
        mjpsi_sqr=3.097*3.097  # GeV^2
        for i in range(len(xlist)):
                x[i]=x[i]+mjpsi_sqr

rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}') 

# Dipole models
# mode label style color
models=[["ipsat", "IPsat",1,1],
	[ "ipsat_nonsatp", "IPsat, nonsatp",2,2],
	[ "iim", "IIM",0,0] ]
# Gluon distributions
gdists=[""]
# xval label width
xvals=[[0.0001,r"10^{-4}",2.2], [0.01, r"10^{-2}",1] ]
minx=9.54
maxx=100
miny=0.2
maxy=1

textsize=14
rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)
style=1

fig=figure()
fig.subplots_adjust(bottom=0.13)
p1=fig.add_subplot(111)
xlabel(r"$M_{J/\Psi}^2\hspace{-0.1} + Q^2$ $[\mathrm{GeV}^2\hspace{-0.3}]$",fontsize=textsize+4)
ylabel(r"$(\mathrm{d} \sigma^A\hspace{-0.4}/\mathrm{d}t)$ $/$ $(A\mathrm{d}\sigma^p\hspace{-0.40}/\mathrm{d}t)$",fontsize=textsize+4)

color=-1
for xbj in xvals: 
    color=color+1
    style=1
    for mode in models:
        x=[]
        y=[]
        filename = mode[0]+"_x" + str(xbj[0]) + ".txt"
        readfile(filename,x,y) 
	AddJpsiMass(x)
        lbl=mode[1] + r", $x_\mathrm{\mathbb{P}}=" + xbj[1] + r"$"
        p1.semilogx(x,y,label=lbl,linestyle=dashes[mode[2]],color=colors[mode[3]], linewidth=xbj[2])
        style=style+1
        if (style>4):
            style=0
  
leg=legend(prop=dict(size=textsize),labelspacing=0.001, columnspacing=-0.03,ncol=2)
leg.draw_frame(False)
fig.suptitle(r"$A=197$, $|t|=0.5$ GeV$^2$", fontsize=textsize+1)
f = "plot_Q.pdf"
print f
axis([minx,maxx,miny,maxy])   
pp = PdfPages(f)
savefig(pp,format='pdf')
pp.close()
