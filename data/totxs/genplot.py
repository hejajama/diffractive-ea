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
	

# Dipole models
models=["ipsat", "ipsat-nofactor", "iim" ]
labels=[r"IPsat", "Non-factorized IPsat", "IIM"]

textsize=16
style=0
minx=9.5
maxx=100
miny=0.3
maxy=80
mjpsi_sqr=3.097*3.097  # GeV^2
style=1
fig=figure()
fig.subplots_adjust(bottom=0.13)
p1=fig.add_subplot(111, xscale="log", yscale="log")
xlabel(r"$M_{J/\Psi}^2 + Q^2 \hspace{0.1} [\mathrm{GeV}^2\hspace{-0.3}]$",fontsize=textsize+4)
ylabel(r"$\sigma_{q\bar q}^p \hspace{0.1} [\mathrm{nb}]$",fontsize=textsize+4)
i=0
col=0
for mode in models:
	x=[]
        y=[]
        filename = mode + ".txt"
        readfile(filename,x,y,1) 
	AddJpsiMass(x)
        lbl=labels[i]
	i=i+1
        p1.plot(x,y,label=lbl,linestyle=dashes[style],color=colors[col],linewidth=2)
        style=style+1
        if (style>4):
            style=0
	col=col+1
	if col>1:
		col=0
    
# Exp. results 
x=[]
y=[]
err=[]
readfile("exp_h1.txt", x, y, 1, err)
AddJpsiMass(x)
lbl="H1"
p1.errorbar(x,y,yerr=err,label=lbl,marker='o',linestyle="none",color=colors[2],linewidth=2)

x=[]
y=[]
err=[]
errorarray=[]
readfile_errorarray("exp_zeus_2err.txt", x, y, 1, errorarray)
AddJpsiMass(x)
lbl="ZEUS"
p1.errorbar(x,y,yerr=errorarray,label=lbl,marker='^',linestyle="none",color=colors[3],linewidth=2)

axis([minx,maxx,miny,maxy])

leg=legend(prop=dict(size=textsize),labelspacing=0.001,loc=1)
leg.draw_frame(False)
fig.suptitle(r"$\gamma^*\hspace{-0.1}p \rightarrow J/\Psi\hspace{0.1}p$, $W=90$ GeV", fontsize=textsize)
print "xs.pdf"
pp = PdfPages("xs.pdf")
savefig(pp,format='pdf')
pp.close()
