#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;
typedef double (*Func) (double,double);

static const int PLOT_STEP = 50;

/* Dirichlet Boundary conditions */
static const double bc_north = 0.0;
static const double bc_south = 0.0; 
static const double bc_west = 4.0;
static const double bc_east = 0.0;

class HeatSolver2D
{
private:
    int nx, ny, nt;
    double tmax, dt;
    double L, dx;
    double H, dy;
    double alfa;
    Func ic;
public:
    HeatSolver2D (int argc, char *argv[]);
    void solve ();
    void computeCoeff (double &a, double &b, double &c);
    void buildMatrix (SpMat &A);
    void setInitialCondition (VectorXd &b);
    void insertCoefficient (int i, int j, double val, vector<T> &coeff);
    void insertBC_Coefficient (vector<T> &coeff);
    void printHeatSolver2D ();
    void writeSolution (int iter, VectorXd x);
    void writeVTK (const vector<VectorXd> sol);
    
    friend ostream& operator<<(ostream&, const HeatSolver2D&);
};

double initialCondition (double x, double y);

void Usage (char *argv[]);