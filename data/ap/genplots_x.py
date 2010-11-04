#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Plot dsigma/dt
import sys
sys.path.append("../")
from matplotlibhelper import *
rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')

# Dipole models
models=[["ipsat", "IPsat"], ["iim", "IIM"] ] #,"proton","iim_BD_5"]
titles=["IPsat", "IIM"]
# Gluon distributions
gdists=[""]
Qvals=[0,10]
#xbj=0.0001

textsize=16
style=1

fig=figure()
p1=fig.add_subplot(111)
xlabel(r"$x_\mathrm{\mathbb{P}}$",fontsize=textsize+4)
ylabel(r"$(\mathrm{d} \sigma^A\hspace{-0.4}/\mathrm{d}t)$ $/$ $(A\mathrm{d}\sigma^p\hspace{-0.4}/\mathrm{d}t)$",fontsize=textsize+4)

i=-1
for q in Qvals:
    i=i+1	    
    for mode in models:
        x=[]
        y=[]
        filename = mode[0]+"_x_q" + str(q) + ".txt"
        readfile(filename,x,y) 
        lbl=mode[1] + r", $Q^2\hspace{-0.4}="+str(q)+r"$ GeV$^2$" 
        p1.semilogx(x,y,label=lbl,linestyle=dashes[style],linewidth=2,color=colors[style])
        style=style+1
        if (style>3):
            style=0
  
   
leg=legend(prop=dict(size=textsize),labelspacing=0.001,loc=2)
leg.draw_frame(False)
fig.suptitle(r"$A=197$, $t=0.5$ GeV$^2$", fontsize=textsize)
f = "plot_x.pdf"
print (f)
    
pp = PdfPages(f)
savefig(pp,format='pdf')
pp.close()
