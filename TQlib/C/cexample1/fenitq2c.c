#include "octqc.h"
extern void c_tqini(int, gtp_equilibrium_data* );
int main(int argc, char **argv) 
{
	gtp_equilibrium_data *ceq;
	// set some defaults
	int n=0;
	char* filename="feni";

	//initilize tq 
	c_tqini(n,ceq);

	return 0;
}

