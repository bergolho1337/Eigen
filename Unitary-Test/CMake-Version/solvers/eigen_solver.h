#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;                              // Triplet --> (lin,col,val)

VectorXd read_rhs (const char rhsName[]);
SpMat read_matrix (const char mName[]);