#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Plot dsigma/dt

from matplotlibhelper import *

# Dipole models
models=["ipsat", "ipnonsat", "iim"]
# Gluon distributions
gdists=["dglap"]
#Q^2 values
Q2vals=[0,3,20]

mint=0
maxt=1.0

textsize=16
style=0

for q in Q2vals:
    fig=figure()
    p1=fig.add_subplot(111)
    xlabel(r"$|t|$ $/$ $\mathrm{GeV}^2$",fontsize=textsize+4)
    ylabel(r"$\mathrm{d} \sigma / \mathrm{d}t$ $/$ $\mathrm{\mu b}/\mathrm{GeV}^2$",fontsize=textsize+4)
    
    for mode in models:
        for gd in gdists:
            x=[]
            y=[]
            filename = mode+"_" + gd + "_" + str(q) + ".txt"
            readfile(filename,x,y,400*1000)  # In units nb/GeV^2

            lbl=mode+ ", " + gd
            p1.semilogy(x,y,label=lbl,marker='-')#,linestyle=dashes[style],linewidth=2)
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
