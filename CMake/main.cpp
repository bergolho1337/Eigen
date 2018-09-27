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
    cout << "Usage:> " << argv[0] << " <matrix_file> <rhs_file> <store>" << endl;
    cout << "<store> = Store way to populate the sparse matrix and the RHS" << endl;
    cout << "\t1 - Using coefficient vector of Triplet(lin,col,val)" << endl;
    cout << "\t2 - Using Unsupported library of Eigen compatible with MatrixMarket format" << endl;
    cout << "=================================================================" << endl;
}

// Convert a matrix from MatrixMarket format to a SparseMatrix<double> type object
SpMat readMatrix (const char mName[])
{
    vector<T> coeff;
    int n, m, k;
    int l, c;
    double val;
    string str;

    ifstream in(mName);
    getline(in,str);
    in >> n >> m >> k;
    for (int i = 0; i < k; i++)
    {
        in >> l >> c >> val;
        l--; c--;
        coeff.push_back(T(l,c,val));
    }
    in.close();

    // Build the sparse matrix by using the coefficient vector
    SpMat A(n,m);
    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();

    return A;
}

// Convert a vector from MatrixMarket format to a VectorXd type object
VectorXd readVector (const char rhsName[])
{
    int n, k;
    string str;
    ifstream in(rhsName);
    getline(in,str);
    in >> n >> k;
    VectorXd b(n);
    for (int i = 0; i < n; i++)
        in >> b(i);
    in.close();

    return b;
}

// MAIN FUNCTION
int main (int argc, char *argv[])
{
    if (argc-1 < 3)
    {
        Usage(argv);
        exit(EXIT_FAILURE);
    }
    else
    {
        SpMat A;
        VectorXd b, x;
        int store = atoi(argv[3]);
        if (store == 1)
        {
            A = readMatrix(argv[1]);
            cout << "\nSparse matrix A" << endl;
            cout << A << endl;

            b = readVector(argv[2]);
            cout << "\nRHS b" << endl;
            cout << b << endl;

        }
        else if (store == 2)
        {
            loadMarket(A,argv[1]);
            cout << "\nSparse matrix A" << endl;
            cout << A << endl;

            loadMarketVector(b,argv[2]);
            cout << "\nRHS b" << endl;
            cout << b << endl;
        }
        
        SparseLU<SpMat> sparseSolver(A);
        x = sparseSolver.solve(b);

        cout << "\nSolution" << endl;
        cout << x << endl;
    }
  
}