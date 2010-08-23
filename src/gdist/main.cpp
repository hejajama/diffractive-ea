#include "GDist.h"
#include <iostream>
#include <sstream>

using namespace std;

int main(int argc, char *argv[])
{
    float xbj, r;
    
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " xbj r" << endl;
        return -1;
    }
    sscanf(argv[1], "%f", &xbj);
    sscanf(argv[2], "%f", &r);
    
    DGLAPDist xg;
    
    cout << "xg(x=" << xbj << ", r="<<r<<") = " << xg.gluedist(xbj,r*r) << endl;
    
    
    return 0;
}

