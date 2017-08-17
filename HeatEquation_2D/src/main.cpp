/*
------------------------------------------------------------------------------------------------------------
    Solve the 2 dimension Heat equation on a quadricular plate.
    - Use Finite Difference Method (Backward in Time Centered in Space - BTCS)
    - Store the matrix related to the linear system on a SparseMatrix structure (Eigen)
    - First decompose the matrix on LU, after that only change the right-hand-side to solve the system.  
------------------------------------------------------------------------------------------------------------   
*/

#include <iostream>
#include "../include/heat.h"
#include "../include/timer.h"

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
        double start, finish, elapsed;
        HeatSolver2D *heat = new HeatSolver2D(argc,argv);
        cout << *heat << endl;

        GET_TIME(start);

        heat->solve();
        
        GET_TIME(finish);
        elapsed = finish - start;
        printf("Elapsed time = %.10lf seconds\n",elapsed);
        delete heat;
    }
    return 0;
}
