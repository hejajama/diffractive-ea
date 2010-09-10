#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Generates a d\sigma/dt plot
# Data is read from dir q{qval}/model.txt


from matplotlibhelper import *

# Dipole models, syntax: [mode:description]
models=[
	["ipsat","IPSat"],
	["ipnonsat","IPNonSat"],
	["iim","IIM"], 
	["ipsat_nonsatp", "IPSat, NonSat p"]	
]

# Q^2
#qsqr=[0,10,50]
qsqr=[0, 10]

textsize=16
style=0

for q in qsqr:
	style=0
	fig=figure()
	p1=fig.add_subplot(111)
	xlabel(r"$|t|$ $/$ $\mathrm{GeV}^2$",fontsize=textsize+4)
	ylabel(r"$\mathrm{d} \sigma^A / \mathrm{d}t$ $/$ $\mathrm{nb}/\mathrm{GeV}^2$",fontsize=textsize+4)
    
	for mode in models:
        	x=[]
	        y=[]
	        filename = "q" + str(q)+"/" + mode[0] + ".txt"
	        if mode=="proton":
	            readfile(filename,x,y,400*1000*197)
	        else:
	            readfile(filename,x,y,400*1000)  # In units nb/GeV^2
	        lbl=mode[1]
	        p1.semilogy(x,y,label=lbl,linestyle=dashes[style],color=colors[style],linewidth=2)
	        style=style+1
	        if (style>4):
	            style=0
     
   
	leg=legend(prop=dict(size=textsize),labelspacing=0.001)
	leg.draw_frame(False)
	fig.suptitle(r"$Q^2=" + str(q) + "$ GeV$^2$, $x=0.0001$", fontsize=textsize)
	f = "q" + str(q) + ".pdf"
	print f    

	pp = PdfPages(f)
	savefig(pp,format='pdf')
	pp.close()

