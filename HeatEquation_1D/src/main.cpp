/*
------------------------------------------------------------------------------------------------------------
    
------------------------------------------------------------------------------------------------------------   
*/

#include <iostream>
#include "../include/heat.h"

using namespace std;

// MAIN FUNCTION
int main (int argc, char *argv[])
{
    if (argc-1 < 5)
    {
        Usage(argv);
        exit(EXIT_FAILURE);
    }
    else
    {
        HeatSolver1D *heat = new HeatSolver1D(argc,argv);
        cout << *heat << endl;

        heat->solve();
    }
    return 0;
}
