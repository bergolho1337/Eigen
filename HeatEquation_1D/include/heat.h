#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;
typedef double (*Func) (double,double);

/* Dirichlet Boundary conditions */
static const double bc_west = 0.0;
static const double bc_east = 0.0;

class HeatSolver1D
{
private:
    int nx;
    double tmax, dt;
    double L, dx;
    double alfa;
    Func ic;
public:
    HeatSolver1D (int argc, char *argv[]);
    void solve ();
    void computeCoeff (double &a, double &b, double &c, double &d);
    void buildMatrix (const double a, const double b, const double c, const double d, SpMat &A);
    void insertCoefficient (int i, int j, double val, vector<T> &coeff);
    void printHeatSolver1D ();
    
    friend ostream& operator<<(ostream&, const HeatSolver1D&);
};

double initialCondition (double x, double L);

void Usage (char *argv[]);