/*
------------------------------------------------------------------------------------------------------------
    Solve the 1 dimension Heat equation on a quadricular plate.
    - Use Finite Difference Method (Implicit)
    - Store the matrix related to the linear system on a SparseMatrix structure (Eigen)
    - First decompose the matrix on LU, after that only change the right-hand-side to solve the system.
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
