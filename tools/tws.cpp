#include "../src/nucleus.h"
#include <iostream>

using namespace std;

int main()
{
	Nucleus nuke(208);
	nuke.Intialize();

	for (double b=0; b<100; b+=0.1)
		cout << b << " " << nuke.T_WS(b) << endl;
	
	return 0;

}
