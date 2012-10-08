#include "pert.h"

int main()
{
	pert p;
	initPert(2, 2, 1, 3, 9, &p);
	//fillPert(&p);
	dumpPdf(&p);
	dumpCdf(&p);
	freePert(&p);
}