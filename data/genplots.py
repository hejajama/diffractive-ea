#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Plot dsigma/dt

from matplotlibhelper import *

# Dipole models
models=["ipsat", "ipsat_nonsatp", "ipnonsat", "iim"]
# Gluon distributions
gdists=["dglap"]
#Q^2 values
Q2vals=[0,10,50]

mint=0
maxt=1.0

textsize=16
style=0

for q in Q2vals:
    style=0
    fig=figure()
    p1=fig.add_subplot(111)
    xlabel(r"$|t|$ $/$ $\mathrm{GeV}^2$",fontsize=textsize+4)
    ylabel(r"$\mathrm{d} \sigma / \mathrm{d}t$ $/$ $\mathrm{nb}/\mathrm{GeV}^2$",fontsize=textsize+4)
    
    for mode in models:
        for gd in gdists:
            x=[]
            y=[]
            filename = "q"+str(q) + "/" + mode + ".txt"
            readfile(filename,x,y,400*1000)  # In units nb/GeV^2

            lbl=mode+ ", " + gd
            p1.semilogy(x,y,label=lbl,linestyle=dashes[style],color=colors[style],linewidth=2)
            style=style+1
            if (style>4):
                style=0
     
    leg=legend(prop=dict(size=textsize),labelspacing=0.001)
    leg.draw_frame(False)
    fig.suptitle(r"$Q^2=\," + str(q) +r"$ GeV$^2$", fontsize=textsize)
    f = "q2_" + str(q) + ".pdf"
    
    pp = PdfPages(f)
    savefig(pp,format='pdf')
    pp.close()
