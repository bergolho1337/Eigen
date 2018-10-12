#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <unsupported/Eigen/SparseExtra>                // Load MatrixMarket format into Eigen
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace Eigen;
using namespace std;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;                              // Triplet --> (lin,col,val)

void Usage (char *argv[])
{
    cout << "=================================================================" << endl;
    cout << "Usage:> " << argv[0] << " <matrix_file> <rhs_file>" << endl;
    cout << "=================================================================" << endl;
}

VectorXd read_rhs (const char rhsName[])
{
    int n, k;
    string str;
    
    ifstream in(rhsName);
    in >> n;
    VectorXd b(n);
    for (int i = 0; i < n; i++)
        in >> b(i);
    in.close();

    return b;
}

SpMat read_matrix (const char mName[])
{
    vector<T> coeff;
    int n, m;
    int l, c;
    double val;
    string str;

    ifstream in(mName);
    in >> n >> m;
    while (in >> l >> c >> val)
    {
        coeff.push_back(T(l,c,val));
    }
    in.close();

    // Build the sparse matrix by using the coefficient vector
    SpMat A(n,m);
    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();

    return A;
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
