/*
 * Simple program which uses DGLAP gluon distribution glass DGLAPDist
 * to calculate Pi^2/2*Nc * alphas(mu(r)^2) * xg(x,mu(r)^2)
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include "gdist_dglap.h"
#include <iostream>
#include <sstream>

using namespace std;

int main(int argc, char *argv[])
{
    float xbj, r;
    
    DGLAPDist xg;
	cout << "#x    xg(x=0.01,q)" << endl;
	for (double x=1; x<1000; x*=1.05)
	{
		cout << x << " " << xg.Gluedist(0.01, 4.0/x) << endl;
	} return 0;

    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " xbj r" << endl;
        return -1;
    }
    sscanf(argv[1], "%f", &xbj);
    sscanf(argv[2], "%f", &r);
    
    
    cout << "xg(x=" << xbj << ", r="<<r<<") = " << xg.Gluedist(xbj,r*r) << endl;
    
    
    return 0;
}

