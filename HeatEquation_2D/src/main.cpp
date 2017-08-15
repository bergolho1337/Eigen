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
    if (argc-1 < 7)
    {
        Usage(argv);
        exit(EXIT_FAILURE);
    }
    else
    {
        HeatSolver2D *heat = new HeatSolver2D(argc,argv);
        cout << *heat << endl;

        heat->solve();
    }
    return 0;
}
