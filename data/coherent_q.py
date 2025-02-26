#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Generates a d\sigma/dt plot
# Data is read from dir q{qval}/model.txt


from matplotlibhelper import *
rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')
textsize=16
rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)

minx=0
maxx=0.26
miny=10
maxy=5e7

# Dipole models, syntax: [mode:description:style:color,width]
models=[
	["ipnonsat","Quasielastic IPnonsat",1,1,2],
	["ipnonsat_coherent", "Coherent IPnonsat",1,1,1],
	["ipsat","IPsat",0,0,2],
        ["ipsat_coherent", "Coherent IPsat",0,0,1],
	["iim", "IIM",3,2,2],
	["iim_coherent", "Coherent IIM",3,2,1]
]

# Q^2
#qsqr=[0,10,50]
qsqr=[0, 10]

style=0

for q in qsqr:
	style=0
	fig=figure()
	p1=fig.add_subplot(111)
	fig.subplots_adjust(bottom=0.11)
	xlabel(r"$|t|$ $[\mathrm{GeV}^2\hspace{-0.3}]$",fontsize=textsize+4)
	ylabel(r"$\mathrm{d} \sigma^A\hspace{-0.4}/\mathrm{d}t$ $[\mathrm{nb}/\mathrm{GeV}^2\hspace{-0.3}]$",fontsize=textsize+4)
    
	for mode in models:
        	x=[]
	        y=[]
	        filename = "q" + str(q)+"/" + mode[0] + ".txt"
	        if mode=="proton":
	            readfile(filename,x,y)
	        else:
	            readfile(filename,x,y)  # In units nb/GeV^2
	        lbl=mode[1]
	        p1.semilogy(x,y,label=lbl,linestyle=dashes[mode[2]],color=colors[mode[3]],linewidth=mode[4])
	        style=style+1
	        if (style>3):
	            style=0
     
   
	axis([minx,maxx,miny,maxy])
	leg=legend(prop=dict(size=textsize),labelspacing=0.001)
	leg.draw_frame(False)
	fig.suptitle(r"$A=197$, $Q^2\hspace{-0.45}=" + str(q) + "$ GeV$^2\hspace{-0.35}$, $x_\mathrm{\mathbb{P}}=0.001$", fontsize=textsize+1)
	f = "coherent_q" + str(q) + ".pdf"
	print f    
	pp = PdfPages(f)
	savefig(pp,format='pdf')
	pp.close()

