#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <unsupported/Eigen/SparseExtra>                // Load MatrixMarket format into Eigen
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;                              // Triplet --> (lin,col,val)

void Usage (const char pName[]);
SpMat readMatrix (const char mName[]);
VectorXd readVector (const char rhsName[]);
void solveProblem (int argc, char *argv[]);