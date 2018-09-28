/*
------------------------------------------------------------------------------------------------------------
    - Solve a sparse linear system given by two files representing the sparse matrix (matrix.mtx) and the
   right-hand-side (rhs.mtx).
    - The matrix and the RHS are provided using the MatrixMarket format.
      (More info about this format: http://math.nist.gov/MatrixMarket/formats.html)
------------------------------------------------------------------------------------------------------------   
*/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include "solver/solver.h"

// MAIN FUNCTION
int main (int argc, char *argv[])
{
    if (argc-1 < 3)
    {
        Usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    else
    {
        cout << "Solving problem ..." << endl;
        solveProblem(argc,argv);
    }
  
}