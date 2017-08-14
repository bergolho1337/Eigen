#include "../include/heat.h"

double initialCondition (double x, double L)
{
    return sin(M_PI*x/L);
}

HeatSolver1D::HeatSolver1D (int argc, char *argv[])
{
    alfa = atof(argv[1]);
    L = atof(argv[2]);
    dx = atof(argv[3]);
    tmax = atof(argv[4]);
    dt = atof(argv[5]);
    ic = initialCondition;
    nx = nearbyint(L/dx);
}

void HeatSolver1D::solve ()
{
    double a, b, c, d;
    computeCoeff(a,b,c,d);
    SpMat A(nx+1,nx+1);
    buildMatrix(a,b,c,d,A);
    cout << A << endl;
}

void HeatSolver1D::computeCoeff (double &a, double &b, double &c, double &d)
{
    a = (1.0/dt)+(2*alfa/(dx*dx));
    b = -alfa/(dx*dx);
    c = b;
    d = 1.0/dt;
}

void HeatSolver1D::insertCoefficient (int i, int j, double val, vector<T> &coeff)
{
    // Make sure that the column index is not out of eange
    if (j >= 0 && j <= nx)
        coeff.push_back(T(i,j,val));
}

void HeatSolver1D::buildMatrix (const double a, const double b, const double c, const double d, SpMat &A)
{
    int n = nx + 1;
    vector<T> coeff;
    insertCoefficient(0,0,1,coeff);             // Boundary coefficient
    for (int i = 1; i < n-1; i++)
    {
        insertCoefficient(i,i,a,coeff);
        insertCoefficient(i,i-1,b,coeff);
        insertCoefficient(i,i+1,c,coeff);
    }
    insertCoefficient(n-1,n-1,1,coeff);             // Boundary coefficient
    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();
}

void HeatSolver1D::printHeatSolver1D ()
{
    cout << "Alfa = " << alfa << endl;
    cout << "L = " << L << endl;
    cout << "dx = " << dx << endl;
    cout << "tmax = " << tmax << endl;
    cout << "dt = " << dt << endl;
}

ostream& operator<<(ostream& ost, const HeatSolver1D& h)
{
    ost << "nx = " << h.nx << endl;
    ost << "Alfa = " << h.alfa << endl;
    ost << "L = " << h.L << endl;
    ost << "dx = " << h.dx << endl;
    ost << "tmax = " << h.tmax << endl;
    ost << "dt = " << h.dt << endl;
    ost << "bc_west = " << bc_west << endl;
    ost << "bc_east = " << bc_east << endl;
    return ost;
}

void Usage (char *argv[])
{
    cout << "=================================================================" << endl;
    cout << "Usage:> " << argv[0] << " <ALPHA> <L> <dx> <tmax> <dt>" << endl;
    cout << "<ALPHA> = Diffusion coefficient" << endl;
    cout << "<L> = Size of the bar" << endl;
    cout << "<dx> = Size of the space discretization" << endl;
    cout << "<tmax> = Maximum simulation time" << endl;
    cout << "<dt> = Size of the time discretization" << endl;
    cout << "=================================================================" << endl;
}