#include <iostream>
#include "solvers/eigen_solver.h"

void Usage (char *argv[])
{
    cout << "=================================================================" << endl;
    cout << "Usage:> " << argv[0] << " <matrix_file> <rhs_file>" << endl;
    cout << "=================================================================" << endl;
}

// MAIN FUNCTION
int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        Usage(argv);
        exit(EXIT_FAILURE);
    }
    else
    {
        SpMat A;
        VectorXd b, x;
        
        // _____________________________
        // Build the matrix
        A = read_matrix(argv[1]);

        // _____________________________
        // Build the right-hand-side
        b = read_rhs(argv[2]);
    
        // _____________________________
        // Build the SparseLU solver
        SparseLU<SpMat> sparseSolver(A);

        // _____________________________
        // Solve the linear system
        x = sparseSolver.solve(b);

        // _____________________________
        // Write the solution to a file
        FILE *out_file = fopen("solution.txt","w+");
        for (int i = 0; i < x.size(); i++)
            fprintf(out_file,"%.10lf\n",x(i));
        fclose(out_file);
    }
  
    return 0;
}
